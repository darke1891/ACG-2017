#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <Eigen/Geometry>
#include <vector>

NORI_NAMESPACE_BEGIN

class PathIntegrator : public Integrator {
public:
    PathIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray, bool is_diffuse, float pb) const {
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Color3f res(0.0f);
        float pl;

        const BSDF *bsdf;
        const Emitter *emitter;
        bsdf = its.mesh->getBSDF();
        emitter = its.mesh->getEmitter();
        BSDFQueryRecord bRec(Vector3f(0.0f, 0.0f, 0.0f));

        const std::vector<Mesh *> meshs = scene->getMeshes();
        const Emitter *emitter2;
        std::vector<int> mesh_idx;
        mesh_idx.clear();
        DiscretePDF dpdf;
        dpdf.clear();
        for (uint32_t idx = 0; idx < meshs.size(); idx++) {
            emitter2 = meshs[idx]->getEmitter();
            if (emitter2 == nullptr)
                continue;
            dpdf.append(1.0f);
            mesh_idx.push_back(idx);
        }
        float light_sum = dpdf.normalize();

        if (emitter != nullptr) {
            Color3f direct_res = emitter->hit(its.p);
            Vector3f d_direct = -ray.d;
            d_direct = its.shFrame.toLocal(d_direct);
            if (d_direct.z() <= 0.0f)
                direct_res *= 0.0f;
            else if (is_diffuse) {
                pl = emitter->get_pdf();
                pl /= light_sum;
                pl /= d_direct.z();
                pl *= its.t * its.t;
                if (pb + pl > 0.0f)
                    direct_res *= pb / (pb + pl);
            }
            res += direct_res;
        }
        else {
            float emitter_u = sampler->next1D();
            float emitter_pdf;
            size_t emitter_id = dpdf.sample(emitter_u, emitter_pdf);
            emitter = meshs[mesh_idx[emitter_id]]->getEmitter();

            Point2f sample = sampler->next2D();
            EmitterSample emitter_sample = emitter->sample(sample);
            Vector3f d = emitter_sample.point - its.p;
            Vector3f d_norm;
            float dis = d.norm();
            d_norm = d;
            d_norm.normalize();
            Ray3f newRay(its.p, d_norm, Epsilon, dis - Epsilon);
            float g_weight = 1.0f / emitter_sample.probability_density;
            if (scene->rayIntersect(newRay))
                g_weight = 0;
            if (emitter_sample.normal.dot(d_norm) >= 0.0f)
                g_weight = 0;
            else
                g_weight *= fabs(its.shFrame.n.dot(d_norm)) * fabs(emitter_sample.normal.dot(d_norm)) / (dis * dis);
            bRec = BSDFQueryRecord(its.shFrame.toLocal(-ray.d), its.shFrame.toLocal(d_norm), ESolidAngle);
            
            Color3f d_color = bsdf->eval(bRec);
            d_color *= g_weight;
            d_color *= emitter_sample.radiance;
            d_color /= emitter_pdf;

            pl = - emitter_sample.normal.dot(d_norm);
            if (pl > 0.0f) {
                pl = emitter->get_pdf() / pl;
                pl *= emitter_pdf;
                pl *= dis * dis;
                pb = bsdf->pdf(bRec);
                d_color *= pl / (pb + pl);
            }
            res += d_color;
        }

        float re_xi = sampler->next1D();
        float xi_threshold = 0.95f;
        if (re_xi < xi_threshold) {
            Point2f re_sample = sampler->next2D();
            bRec = BSDFQueryRecord(its.shFrame.toLocal(-ray.d));
            Color3f re_color = bsdf->sample(bRec, re_sample);
            if ((re_color.r() != 0.0f) || (re_color.g() != 0.0f) || (re_color.b() != 0.0f)) {
                pb = bsdf->pdf(bRec);
                Ray3f new_ray(its.p, its.shFrame.toWorld(bRec.wo));
                Color3f re_Lix = Li(scene, sampler, new_ray, bsdf->isDiffuse(), pb);
                res += re_Lix * re_color / xi_threshold * bRec.eta * bRec.eta;
            }
        }

        return res;
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        return Li(scene, sampler, ray, false, 0.0f);
    }

    std::string toString() const {
        return "PathIntegrator[]";
    }

};

NORI_REGISTER_CLASS(PathIntegrator, "path");
NORI_NAMESPACE_END