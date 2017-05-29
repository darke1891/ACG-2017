#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/dpdf.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

struct VolumeHeteroFlameData {
    Transform toWorld, toLocal;
    float*** flame = nullptr;
    int* flame_index = nullptr;
    float*** heat = nullptr;
    int* heat_index = nullptr;
};

class VolumeEmitterHetero : public Emitter {
public:
    VolumeEmitterHetero(const PropertyList &props) {
        radiance = props.getColor("radiance");
    }

    void transfer_data(VolumeHeteroFlameData &flameData);

    std::string toString() const {
        return "VolumeEmitterHetero";
    }
    virtual EmitterSample sample(Point2f &p) const;
    virtual EmitterSample sample(Point2f &p, float p2) const;
    virtual Color3f hit(Point3f p) const;
    virtual float get_pdf(Point3f p) const;

private:
    Color3f radiance;
    Transform toWorld, toLocal;
    float*** flame = nullptr;
    int* flame_index = nullptr;
    float*** heat = nullptr;
    int* heat_index = nullptr;
};

NORI_NAMESPACE_END