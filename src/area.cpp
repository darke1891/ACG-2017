#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/dpdf.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &props) {
        radiance = props.getColor("radiance");
    }

    std::string toString() const {
        return "AreaLightEmitter";
    }
    virtual void emit_activate(Mesh *);
    virtual EmitterSample sample(Point2f &p) const;
    virtual Color3f hit(Point3f p) const;
    virtual float get_pdf(Point3f p) const;

private:
    Color3f radiance;
    Mesh *parent_mesh = nullptr;
    DiscretePDF dpdf;
    float probability_density;
};

Color3f AreaLight::hit(Point3f p) const {
    return radiance;
}

float AreaLight::get_pdf(Point3f p) const {
    return probability_density;
}

void AreaLight::emit_activate(Mesh *parent) {
    parent_mesh = parent;
    dpdf.clear();
    for (uint32_t idx = 0; idx < parent_mesh->getTriangleCount(); ++idx)
        dpdf.append(parent_mesh->surfaceArea(idx));
    probability_density = dpdf.normalize();
    probability_density = 1.0f / probability_density;
}

EmitterSample AreaLight::sample(Point2f &p) const {
    float u, v;
    u = p.x();
    v = p.y();
    size_t idx = dpdf.sampleReuse(u);
    float a, b;
    a = 1 - sqrt(1 - u);
    b = v * sqrt(1 - u);
    const MatrixXu &m_F = parent_mesh->getIndices();
    const MatrixXf &m_V = parent_mesh->getVertexPositions();
    const MatrixXf &m_N = parent_mesh->getVertexNormals();
    uint32_t i0 = m_F(0, idx), i1 = m_F(1, idx), i2 = m_F(2, idx);
    Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);
    Point3f point = p0 + (p1 - p0) * a + (p2 - p0) * b;
    Vector3f normal;
    if (m_N.size() <= 0) {
        Vector3f e0 = p1 - p0;
        Vector3f e1 = p2 - p1;
        normal = e0.cross(e1);
    }
    else
        normal = m_N.col(i0) * (1 - a - b) + m_N.col(i1) * a + m_N.col(i2) * b;
    normal = normal.normalized();

    EmitterSample res;
    res.probability_density = probability_density;
    res.point = point;
    res.normal = normal;
    res.radiance = radiance;
    res.is_surface = true;
    return res;
}

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END