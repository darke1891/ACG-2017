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

Scene::Scene(const PropertyList &props) {
    m_accel = new Accel();

    volumeMedia.clear();
    m_vacuummedia = (VolumeMedia *)(NoriObjectFactory::createInstance("volumemediavacuum", props));
    volumeMedia.insert(std::pair<std::string, const VolumeMedia*>("", m_vacuummedia));
}

Scene::~Scene() {
    delete m_accel;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;
    for (auto it = volumeMedia.begin(); it != volumeMedia.end(); ++it)
        delete (it->second);
    volumeMedia.clear();
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

    for (auto t_mesh: m_meshes) {
        VolumeSurface *t_volumesurface = t_mesh->getNCVolumeSurface();
        if (t_volumesurface)
            t_volumesurface->activate(volumeMedia);
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

        case EVolumeMedia: {
            VolumeMedia* vMedia = static_cast<VolumeMedia *>(obj);
            std::string t_name = vMedia->getName();
            if (volumeMedia.find(t_name) == volumeMedia.end()) {
                volumeMedia.insert(std::pair<std::string, const VolumeMedia*>(t_name, vMedia));
            }
            else
                throw NoriException("Repeated name of volume media!");
            break;
        }
        
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
    float dis_total = 0.0f;
    Ray3f m_ray(ray);

    bool res_hit = false;

    while (true) {
        res_hit = m_accel->rayIntersect(m_ray, its, false);
        if ((res_hit) && (its.mesh)) {
            const VolumeSurface* vSurface = its.mesh->getVolumeSurface();
            if (vSurface) {
                const VolumeMedia* vMedia = vSurface->dirMedia(its.shFrame.toLocal(ray.d));
                float dis = vMedia->sample_dis(m_ray, its.t, sampler);
                if (dis < its.t) {
                    its = Intersection();
                    its.t = dis;
                    its.p = m_ray.o + dis * m_ray.d;
                    its.shFrame = Frame((-ray.d).normalized());
                    its.geoFrame = its.shFrame;
                    its.m_media = vMedia;
                    its.is_surface = false;
                    break;
                }
            }
            if (its.mesh->getBSDF())
                break;
            else if (vSurface) {
                dis_total += its.t;
                m_ray.maxt -= its.t;
                m_ray.o += m_ray.d * its.t;
            }
            else
                break;
        }
        else
            break;
    }
    its.t += dis_total;
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
    return res_hit;
}

bool Scene::rayIntersect(const Ray3f &ray, Sampler *sampler) const {
    Intersection its;
    return rayIntersect(ray, its, sampler);
}

const BoundingBox3f &Scene::getBoundingBox() const {
    return m_accel->getBoundingBox();
}

std::vector<const Emitter*> Scene::get_emitterMedia() const {
    std::vector<const Emitter*> res;
    res.clear();
    for (auto it = volumeMedia.begin(); it != volumeMedia.end(); ++it) {
        const Emitter* emitter = it->second->getEmitter();
        if (emitter)
            res.push_back(emitter);
    }
    return res;
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
