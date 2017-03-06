/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <nori/timer.h>
#include "tbb/parallel_invoke.h"
#include <algorithm>

NORI_NAMESPACE_BEGIN

#define ACCEL_OCTREE_ON 1
#define ACCEL_OCTREENODE_IMPROVED_TRAVERSAL 1
#define ACCEL_OCTREE_PARALLEL_CONSTRUCTION 0

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() {
    /* Nothing to do here for now */
    Timer timer;
    m_octree = new OctreeNode();
    m_octree->setMesh(m_mesh);
    #if ACCEL_OCTREE_PARALLEL_CONSTRUCTION == 1
        m_octree->splitNode(m_mesh, 0, false);
    #else
        m_octree->splitNode(m_mesh, 0);
    #endif

    cout << "It took " << timer.elapsedString() << " to construct octree" << endl;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)



    #if ACCEL_OCTREE_ON == 0
        /* Brute force search through all triangles */
        for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
            float u, v, t;
            if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
                /* An intersection was found! Can terminate
                   immediately if this is a shadow ray query */
                if (shadowRay)
                    return true;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_mesh;
                f = idx;
                foundIntersection = true;
            }
        }
    #else
        foundIntersection = m_octree->rayIntersect(m_mesh, ray, its, f, shadowRay);
        if (shadowRay && foundIntersection)
            return true;
    #endif

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

bool Accel::OctreeNode::rayIntersect(Mesh *mesh, Ray3f &ray, Intersection &its, uint32_t &f, bool shadowRay) const{
    bool foundIntersection = false;  // Was an intersection found so far?

    if (children != nullptr) {
        std::pair<float, uint32_t> to_search[8];
        float u, v;
        for (uint32_t idx = 0; idx < 8; idx++) {
            to_search[idx] = std::make_pair(ray.maxt, 8);
            if (children[idx].m_bbox.rayIntersect(ray, u, v))
                if (ray.mint <= v && u <= ray.maxt)
                    to_search[idx] = std::make_pair(u, idx);
        }

        #if ACCEL_OCTREENODE_IMPROVED_TRAVERSAL == 1
            std::sort(to_search, to_search + 8);
        #endif

        for (uint32_t i = 0; i < 8; i++)
            if (to_search[i].second < 8)
                if (children[to_search[i].second].rayIntersect(mesh, ray, its, f, shadowRay)) {
                    #if ACCEL_OCTREENODE_IMPROVED_TRAVERSAL == 1
                        return true;
                    #else
                        if (shadowRay)
                            return true;
                        foundIntersection = true;
                    #endif
                }
    }

    if (triangles != nullptr)
        for (uint32_t i = 0; i < num_triangle; i++) {
            uint32_t idx = triangles[i];
            float u, v, t;
            if (mesh->rayIntersect(idx, ray, u, v, t)) {
                /* An intersection was found! Can terminate
                   immediately if this is a shadow ray query */
                if (shadowRay)
                    return true;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = mesh;
                f = idx;
                foundIntersection = true;
            }
        }

    return foundIntersection;
}

void Accel::OctreeNode::splitNode (Mesh *mesh, uint32_t depth, bool full) {
    if (temp_triangles == nullptr)
        return;

    uint32_t num = temp_triangles->size();
    if ((num < 10) || (depth >= 10)) {
        makeLeaf();
        return;
    }
    
    children = new OctreeNode[8];
    bool sameChild = true;

    for (uint32_t i=0; i<8; i++) {
        BoundingBox3f child_bbox(m_bbox.getCorner(i));
        child_bbox.expandBy(m_bbox.getCenter());
        std::vector<uint32_t> *child_triangles = new std::vector<uint32_t>;
        child_triangles->clear();
        for (auto idx: *temp_triangles)
            if (child_bbox.overlaps(mesh->getBoundingBox(idx)))
                child_triangles->push_back(idx);
        if (child_triangles->size() != temp_triangles->size())
            sameChild = false;
        children[i].m_bbox = child_bbox;
        children[i].temp_triangles = child_triangles;
    }
    if (sameChild) {
        delete []children;
        children = nullptr;
        makeLeaf();
    }
    else {
        nullTemp();
        if (full || num < 50)
            for (uint32_t i=0; i<8; i++)
                children[i].splitNode(mesh, depth + 1, true);
        else {
            tbb::parallel_invoke([=]{children[0].splitNode(mesh, depth + 1, full);},
                                 [=]{children[1].splitNode(mesh, depth + 1, full);},
                                 [=]{children[2].splitNode(mesh, depth + 1, full);},
                                 [=]{children[3].splitNode(mesh, depth + 1, full);},
                                 [=]{children[4].splitNode(mesh, depth + 1, full);},
                                 [=]{children[5].splitNode(mesh, depth + 1, full);},
                                 [=]{children[6].splitNode(mesh, depth + 1, full);},
                                 [=]{children[7].splitNode(mesh, depth + 1, full);}
                                 );
        }
    }
}

void Accel::OctreeNode::makeLeaf() {
    num_triangle = temp_triangles->size();
    if (num_triangle > 0) {
        triangles = new uint32_t[num_triangle];
        for (uint32_t i=0; i<num_triangle; i++)
            triangles[i] = (*temp_triangles)[i];
    }
    nullTemp();
}

void Accel::OctreeNode::setMesh (Mesh *mesh) {
    temp_triangles = new std::vector<uint32_t>;
    temp_triangles->clear();
    for (uint32_t idx = 0; idx < mesh->getTriangleCount(); ++idx)
        temp_triangles->push_back(idx);
    m_bbox = mesh->getBoundingBox();
}


void Accel::OctreeNode::nullTemp() {
    delete temp_triangles;
    temp_triangles = nullptr;
}


NORI_NAMESPACE_END

