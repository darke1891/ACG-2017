#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/volumemedia.h>
#include <Eigen/Geometry>
#include <vector>

NORI_NAMESPACE_BEGIN

class PathLoopIntegrator : public Integrator {
public:
    PathLoopIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {

        Ray3f m_ray(ray);
        bool is_diffuse = false;
        Color3f res(0.0f);
        Color3f now_one(1.0f);
        float re_xi = 0.0f;
        float xi_threshold_init = 0.95f;
        float xi_threshold = xi_threshold_init;
        int jump = 0;

        const std::vector<Mesh *> meshs = scene->getMeshes();
        std::vector<const Emitter *> emitters;
        emitters.clear();
        const Emitter *emitter2;
        DiscretePDF dpdf;
        dpdf.clear();
        for (uint32_t idx = 0; idx < meshs.size(); idx++) {
            emitter2 = meshs[idx]->getEmitter();
            if (emitter2 == nullptr)
                continue;
            dpdf.append(1.0f);
            emitters.push_back(emitter2);
        }
        const SceneBox *sbox = scene->get_scenebox();
        if (sbox)
            if (sbox->getEmitter()) {
                dpdf.append(1.0f);
                emitters.push_back(sbox->getEmitter());
            }
        dpdf.normalize();

        while (true) {
            Intersection its;
            if (!scene->rayIntersect(m_ray, its, sampler))
                break;

            xi_threshold = std::min(xi_threshold_init, now_one.maxValue() * 10.0f);

            const BSDF *bsdf;
            const Emitter *emitter;
            bsdf = its.getBSDF();
            emitter = its.getEmitter();
            BSDFQueryRecord bRec(Vector3f(0.0f, 0.0f, 0.0f), sampler);

            if (!is_diffuse) {
                if (emitter != nullptr) {
                    Color3f direct_res = emitter->hit(its.p);
                    Vector3f d_direct = -m_ray.d;
                    d_direct = its.shFrame.toLocal(d_direct);
                    if (d_direct.z() <= 0.0f)
                        direct_res *= 0.0f;
                    res += direct_res * now_one;
                }
                bRec = BSDFQueryRecord(its.shFrame.toLocal(-m_ray.d), sampler);
            }
            else if ((bsdf) && (bsdf->isDiffuse())) {
                float emitter_u = sampler->next1D();
                float emitter_pdf;
                size_t emitter_id = dpdf.sample(emitter_u, emitter_pdf);
                emitter = emitters[emitter_id];

                Point2f sample = sampler->next2D(2);
                EmitterSample emitter_sample = emitter->sample(sample);
                Vector3f d = emitter_sample.point - its.p;
                Vector3f d_norm;
                float dis = d.norm();
                float short_dis = 0.01f;
                if (dis < short_dis)
                    dis = short_dis;
                d_norm = d;
                d_norm.normalize();
                Ray3f newRay(its.p, d_norm, Epsilon, dis - Epsilon);
                float g_weight = 1.0f / emitter_sample.probability_density;

                bool zero_light = false;

                if (scene->rayIntersect(newRay, sampler))
                    zero_light = true;

                if (emitter_sample.normal.dot(d_norm) >= 0.0f)
                    zero_light = true;
                else
                    g_weight *= fabs(emitter_sample.normal.dot(d_norm)) / (dis * dis);

                if (its.is_surface){
                    g_weight *= fabs(its.shFrame.n.dot(d_norm));
                }
                else
                    if (dis < 0.0f) {
                        cout << dis << endl;
                        cout << emitter_sample.point << endl;
                        cout << its.p << endl;
                        cout << "---------" << endl;
                    }

                bRec = BSDFQueryRecord(its.shFrame.toLocal(-m_ray.d), its.shFrame.toLocal(d_norm), ESolidAngle, sampler);
                if (!zero_light) {
                    Color3f d_color = bsdf->eval(bRec);
                    d_color *= g_weight;
                    d_color *= emitter_sample.radiance;
                    d_color /= emitter_pdf;
                    res += d_color * now_one;
                }
            }

            re_xi = sampler->next1D();
            if ((re_xi < xi_threshold) || (jump < 4)) {
                if (bsdf) {
                    Point2f re_sample = sampler->next2D(3 + jump);
                    Color3f re_color = bsdf->sample(bRec, re_sample);
                    if (re_color.maxValue() > Epsilon) {
                        m_ray = Ray3f(its.p, its.shFrame.toWorld(bRec.wo));
                        is_diffuse = bsdf->isDiffuse();
                        if (is_diffuse)
                            jump++;
                        now_one *= re_color / xi_threshold;
                    }
                    else
                        break;
                }
                else {
                    m_ray = Ray3f(its.p, its.shFrame.toWorld(- bRec.wi));
                    is_diffuse = false;
                    now_one /= xi_threshold;
                }
            }
            else
                break;
        }

        return res;
    }

    std::string toString() const {
        return "PathLoopIntegrator[]";
    }

};

NORI_REGISTER_CLASS(PathLoopIntegrator, "loop");
NORI_NAMESPACE_END