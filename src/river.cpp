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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class River : public BSDF {
public:
    River(const PropertyList &props) { 
        m_albedo = props.getColor("albedo", Color3f(0.5f));
        reflect= props.getFloat("reflect", 0.0f);
    }

    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        return m_albedo * INV_PI * (1.0f - reflect);
    }

    float pdf(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;
        return INV_PI * fabs(Frame::cosTheta(bRec.wo)) * (1.0f - reflect);
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        float px, py;
        px = sample.x();
        py = sample.y();

        bRec.measure = ESolidAngle;

        if (px < reflect) {
            bRec.wo = Vector3f(
                -bRec.wi.x(),
                -bRec.wi.y(),
                 bRec.wi.z()
            );
            bRec.measure = EDiscrete;

            /* Relative index of refraction: no change */
            bRec.eta = 1.0f;

            return Color3f(1.0f);
        }

        px = (px - reflect) / (1.0f - reflect);

        bRec.wo = Warp::squareToCosineHemisphere(Point2f(px, py));

        if (Frame::cosTheta(bRec.wi) <= 0)
            bRec.wo *= -1;

        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        /* eval() / pdf() * cos(theta) = albedo. There
           is no need to call these functions. */
        return m_albedo;
    }

    bool isDiffuse() const {
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Bidiffuse[\n"
            "  albedo = %s\n"
            "  reflect= %f\n"
            "]", m_albedo.toString(),
            reflect);
    }
private:
    float reflect;
    Color3f m_albedo;
};

NORI_REGISTER_CLASS(River, "river");
NORI_NAMESPACE_END
