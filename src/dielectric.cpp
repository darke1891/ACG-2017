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
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        bRec.measure = EDiscrete;

        if (Frame::cosTheta(bRec.wi) == 0.0f)
            return Color3f(0.0f);

        float eta1, eta2;
        Vector3f n;
        if (Frame::cosTheta(bRec.wi) <= 0.0f) {
            eta1 = m_intIOR;
            eta2 = m_extIOR;
            n = Vector3f(0, 0, -1.0f);
        }
        else {
            eta1 = m_extIOR;
            eta2 = m_intIOR;
            n = Vector3f(0, 0, 1.0f);
        }
        float F = fresnel(fabs(Frame::cosTheta(bRec.wi)), eta1, eta2);
        if (sample.x() > F) {
            float eta = eta1 / eta2;
            float weight0 = bRec.wi.dot(n);
            weight0 = sqrt(1.0f - eta * eta * (1 - weight0 * weight0));
            Vector3f wt = - weight0 * n;
            wt -= eta * (bRec.wi - bRec.wi.dot(n) * n);

            bRec.wo = wt;
            bRec.eta = eta;
        }
        else {
            bRec.wo = 2 * n * bRec.wi.dot(n) - bRec.wi;
            bRec.eta = 1.0f;
        }
        if (Frame::cosTheta(bRec.wo) == 0.0f)
            return Color3f(0.0f);

        return Color3f(bRec.eta * bRec.eta);
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
