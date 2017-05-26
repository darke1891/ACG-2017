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

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/emitter.h>
#include <nori/scenebox.h>

NORI_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &) {
    m_accel = new Accel();

    PropertyList m_list;
    m_bsdf = (BSDF *)(NoriObjectFactory::createInstance("volumnHG", m_list));
    theta_t = 1.0f;
}

Scene::~Scene() {
    delete m_accel;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;
}

void Scene::activate() {
    m_accel->build();

    if (!m_integrator)
        throw NoriException("No integrator was specified!");
    if (!m_camera)
        throw NoriException("No camera was specified!");
    
    if (!m_sampler) {
        /* Create a default (independent) sampler */
        m_sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
    }

    cout << endl;
    cout << "Configuration: " << toString() << endl;
    cout << endl;
}

void Scene::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EMesh: {
                Mesh *mesh = static_cast<Mesh *>(obj);
                m_accel->addMesh(mesh);
                m_meshes.push_back(mesh);
            }
            break;
        
        case EEmitter: {
                // Emitter *emitter = static_cast<Emitter *>(obj);
                /* TBD */
                throw NoriException("Scene::addChild(): You need to implement this for emitters");
            }
            break;
        case ESceneBox: {
            if (m_box)
                throw NoriException("There can only be one scene box!");
            m_box = static_cast<SceneBox *>(obj);
            break;
        }

        case ESampler:
            if (m_sampler)
                throw NoriException("There can only be one sampler per scene!");
            m_sampler = static_cast<Sampler *>(obj);
            break;

        case ECamera:
            if (m_camera)
                throw NoriException("There can only be one camera per scene!");
            m_camera = static_cast<Camera *>(obj);
            break;
        
        case EIntegrator:
            if (m_integrator)
                throw NoriException("There can only be one integrator per scene!");
            m_integrator = static_cast<Integrator *>(obj);
            break;

        default:
            throw NoriException("Scene::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
    }
}

bool Scene::rayIntersect(const Ray3f &ray, Intersection &its, Sampler *sampler) const {
    bool res_hit = m_accel->rayIntersect(ray, its, false);
    if (m_box) {
        Intersection its_box;
        bool box_hit = m_box->rayIntersect(ray, its_box);
        if (box_hit) {
            if ((!res_hit) || (its_box.t < its.t)) {
                res_hit = box_hit;
                its = its_box;
            }
        }
    }
    if ((m_bsdf) && (theta_t > 0.0f)) {
        float dis = sampler->next1D();
        dis = -(log(dis) / theta_t);
        if ((res_hit) && ((dis < its.t) && (dis < ray.maxt))) {
            res_hit = true;
            its = Intersection();
            its.t = dis;
            its.p = ray.o + dis * ray.d;
            its.shFrame = Frame((- ray.d).normalized());
            its.geoFrame = its.shFrame;
            its.m_bsdf = m_bsdf;
            its.is_surface = false;
        }
    }
    return res_hit;
}

bool Scene::rayIntersect(const Ray3f &ray, Sampler *sampler) const {
    Intersection its;
    return rayIntersect(ray, its, sampler);
}

const BoundingBox3f &Scene::getBoundingBox() const {
    return m_accel->getBoundingBox();
}

std::string Scene::toString() const {
    std::string meshes;
    for (size_t i=0; i<m_meshes.size(); ++i) {
        meshes += std::string("  ") + indent(m_meshes[i]->toString(), 2);
        if (i + 1 < m_meshes.size())
            meshes += ",";
        meshes += "\n";
    }

    return tfm::format(
        "Scene[\n"
        "  integrator = %s,\n"
        "  sampler = %s\n"
        "  camera = %s,\n"
        "  meshes = {\n"
        "  %s  }\n"
        "]",
        indent(m_integrator->toString()),
        indent(m_sampler->toString()),
        indent(m_camera->toString()),
        indent(meshes, 2)
    );
}

NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END
