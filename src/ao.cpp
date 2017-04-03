#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AoIntegrator : public Integrator {
public:
    AoIntegrator(const PropertyList &props) {
        /* No parameters this time */
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Vector3f n, n1, n2;
        n = its.shFrame.n;
        n1 = n;
        if (fabs(n.x()) > 0.2)
            n1.y() += 1;
        else
            n1.x() += 1;
        n1 = n1.normalized();
        n1 -= n * n.dot(n1);
        n1 = n1.normalized();
        n2.x() = n1.y() * n.z() - n1.z() * n.y();
        n2.y() = n1.z() * n.x() - n1.x() * n.z();
        n2.z() = n1.x() * n.y() - n1.y() * n.x();

        Point2f sample = sampler->next2D();
        Vector3f d = Warp::squareToCosineHemisphere(sample);
        d = d.x() * n1 + d.y() * n2 + d.z() * n;

        Ray3f newRay(its.p, d);
        if (scene->rayIntersect(newRay)) {
            return Color3f(0.0f);
        }
        Color3f res = Color3f(1.0f);
        return res;
    }

    std::string toString() const {
        return "AoIntegrator[]";
    }
};

NORI_REGISTER_CLASS(AoIntegrator, "ao");
NORI_NAMESPACE_END