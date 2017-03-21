#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList &props) {
        position = props.getPoint("position");
        energy = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Vector3f d = position - its.p;
        double dis = d.norm();
        d = d.normalized();
        Ray3f newRay(its.p, d, Epsilon, dis);
        if (scene->rayIntersect(newRay)) {
            return Color3f(0.0f);
        }
        double w = 0;
        double cos_theta = its.shFrame.n.dot(d);
        if (cos_theta > 0)
            w = cos_theta;
        w /= dis * dis * 4 * M_PI * M_PI;
        return energy * w;
    }

    std::string toString() const {
        return "SimpleIntegrator[]";
    }

    Point3f position;
    Color3f energy;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END