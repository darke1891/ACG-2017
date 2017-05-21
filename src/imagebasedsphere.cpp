#include <nori/emitter.h>
#include <nori/scenebox.h>
#include <nori/mesh.h>
#include <nori/dpdf.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/hierarchical.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

class ImageBasedSphereEmitter : public Emitter {
public:
    ImageBasedSphereEmitter(const PropertyList &props);
    ~ImageBasedSphereEmitter();
    std::string toString() const {
        return "ImageBasedSphereEmitter";
    }
    virtual EmitterSample sample(Point2f &p) const;
    virtual Color3f hit(Point3f p) const;
    virtual float get_pdf(Point3f p) const;
    virtual bool for_scene() const;
private:
    Point3f center;
    float radius;
    float scale;
    HierarchicalSampler *hSampler = nullptr;
};

class ImageBasedSphere : public SceneBox {
public:
    ImageBasedSphere(const PropertyList &props);
    ~ImageBasedSphere();
    virtual const Emitter *getEmitter() const;
    virtual const BSDF *getBSDF() const;
    virtual bool rayIntersect(const Ray3f &ray, Intersection &its) const;

    std::string toString() const {
        return "ImageBasedSphere";
    }
    virtual void activate();

private:
    BSDF *m_bsdf = nullptr;
    Point3f center;
    float radius;
    ImageBasedSphereEmitter *m_emitter;
};

ImageBasedSphere::ImageBasedSphere(const PropertyList &props) {
    center = props.getPoint("center");
    radius = props.getFloat("radius");
    std::string file_name = props.getString("light-image", "");
    if (file_name != "") {
        m_emitter = new ImageBasedSphereEmitter(props);
    }
}

ImageBasedSphere::~ImageBasedSphere() {
    delete m_bsdf;
    if (m_emitter)
        delete m_emitter;
}

bool ImageBasedSphere::rayIntersect(const Ray3f &ray, Intersection &its) const {
    Point3f oc = ray.o - center;
    float a, bh, c, t;
    a = ray.d.norm();
    a = a * a;
    bh = ray.d.dot(oc);
    c = oc.norm();
    c = c * c - radius * radius;
    c = bh * bh - c;
    if (c >= 0) {
        c = sqrt(c);
        t = - c - bh;
        if (t > ray.mint)
            return false;
        if (t < ray.mint) {
            t = c - bh;
            if ((t < ray.mint) || (t > ray.maxt))
                return false;
        }
        its.t = t;
        its.p = ray.o + t * ray.d;
        its.shFrame = Frame((center - its.p).normalized());
        its.geoFrame = its.shFrame;
        its.sbox = this;
        return true;
    }
    else
        return false;
}

const Emitter *ImageBasedSphere::getEmitter() const {
    return m_emitter;
}

const BSDF *ImageBasedSphere::getBSDF() const {
    return m_bsdf;
}

void ImageBasedSphere::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }
}

ImageBasedSphereEmitter::ImageBasedSphereEmitter(const PropertyList &props) {
    center = props.getPoint("center");
    radius = props.getFloat("radius");
    scale = props.getFloat("scale", 1.0f);
    std::string file_name = props.getString("light-image", "");
    hSampler = new HierarchicalSampler();
    hSampler->setImage(file_name);
}

ImageBasedSphereEmitter::~ImageBasedSphereEmitter() {
    delete hSampler;
}

EmitterSample ImageBasedSphereEmitter::sample(Point2f &p) const {
    EmitterSample sample;
    Point2f pick = hSampler->squareToHierarchical(p);
    sample.probability_density = hSampler->squareToHierarchicalPdf(pick);
    sample.probability_density /= radius * radius * M_PI * 4.0f;
    sample.normal.y() = pick.y() * 2.0f - 1.0f;
    float xz = sqrt(1.0f - sample.normal.y() * sample.normal.y());
    sample.normal.x() = - cos(pick.x() * 2.0f * M_PI) * xz;
    sample.normal.z() = - sin(pick.x() * 2.0f * M_PI) * xz;
    sample.normal *= -1.0f;
    sample.normal.normalize();

    sample.point = center - sample.normal * radius;
    sample.radiance = hSampler->hit(pick) * scale;
    return sample;
}

Color3f ImageBasedSphereEmitter::hit(Point3f p) const {
    Point3f po = p - center;
    po.normalize();
    Point2f p2;
    p2.y() = (po.y() + 1.0f) / 2.0f;
    po.y() = 0.0f;
    po.normalize();
    p2.x() = acos(- po.x());
    if (po.z() > 0.0f)
        p2.x() = 2.0f * M_PI - p2.x();
    p2.x() /= 2.0f * M_PI;
    return hSampler->hit(p2) * scale;
}

float ImageBasedSphereEmitter::get_pdf(Point3f p) const {
    Point3f po = p - center;
    po.normalize();
    Point2f p2;
    p2.y() = (po.y() + 1.0f) / 2.0f;
    po.y() = 0.0f;
    po.normalize();
    p2.x() = acos(- po.x());
    if (po.z() > 0.0f)
        p2.x() = 2.0f * M_PI - p2.x();
    p2.x() /= 2.0f * M_PI;
    float res = hSampler->squareToHierarchicalPdf(p2);
    res /= radius * radius * M_PI * 4.0f;
    return res;
}

bool ImageBasedSphereEmitter::for_scene() const {
    return true;
}

NORI_REGISTER_CLASS(ImageBasedSphere, "image-based-sphere");
NORI_NAMESPACE_END