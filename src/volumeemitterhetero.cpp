#include <nori/volumeemitterhetero.h>

NORI_NAMESPACE_BEGIN

Color3f VolumeEmitterHetero::hit(Point3f p) const {
    Point3f pos;
    float flame_value, heat_value, density_value;
    pos = toLocal * p;

    density_value = interpolation_3d(density, density_index, pos);
    flame_value = interpolation_3d(flame, flame_index, pos);
    if (flame_value <= 0.0f)
        return Color3f(0.0f);
    heat_value = interpolation_3d(heat, heat_index, pos);
    if (heat_value <= 0.0f)
        return Color3f(0.0f);

    if (density_value > 0)
        flame_value = flame_value / (flame_value + density_value);
    else
        flame_value = 1.0f;

    float r,g,b;
    r = radiance.r() * flame_value;
    g = radiance.g() * flame_value;
    b = radiance.b() * flame_value;
    r *= std::max(std::min(heat_value * 1.25f - 1.1f, 1.0f), 0.0f);
    g *= std::max(std::min(heat_value * 2.5f - 2.8f, 1.0f), 0.0f);
    b *= std::max(std::min(heat_value * 5.0f - 6.8f, 1.0f), 0.0f);
    return Color3f(r, g, b);
}

float VolumeEmitterHetero::get_pdf(Point3f p) const {
    return 1.0f;
}

EmitterSample VolumeEmitterHetero::sample(Point2f &p) const {
    EmitterSample res;
    return res;
}

EmitterSample VolumeEmitterHetero::sample(Point2f &p, float p2) const {
    EmitterSample res;
    Point3f pos = Point3f(p.x(), p.y(), p2);
    res.point = toWorld * pos;
    res.radiance = hit(res.point);
    res.probability_density = get_pdf(res.point);
    res.is_surface = false;
    return res;
}

void VolumeEmitterHetero::transfer_data(VolumeHeteroFlameData &flameData) {
    toWorld = flameData.toWorld;
    toLocal = flameData.toLocal;
    density = flameData.density;
    density_index = flameData.density_index;
    flame = flameData.flame;
    flame_index = flameData.flame_index;
    heat = flameData.heat;
    heat_index = flameData.heat_index;
}


NORI_NAMESPACE_END