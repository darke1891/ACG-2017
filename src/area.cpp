#include <nori/emitter.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
public:
    AreaLight(const PropertyList &props) {
        radiance = props.getColor("radiance");
    }

    std::string toString() const {
        return "AreaLightEmitter";
    }
    void setParent(NoriObject *);
    void activate();
    ~AreaLight();
private:
    Color3f radiance;
    Mesh *m_parent;
};

void AreaLight::setParent(NoriObject *parent) {
    m_parent = (Mesh *)parent;
}

void AreaLight::activate() {

}

AreaLight::~AreaLight() {
    delete m_parent;
}

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END