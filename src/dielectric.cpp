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

        float eta1, eta2;
        Vector3f n;
        if (Frame::cosTheta(bRec.wi) <= 0) {
            eta1 = m_intIOR;
            eta2 = m_extIOR;
            n = Vector3f(0, 0, -1);
        }
        else {
            eta1 = m_extIOR;
            eta2 = m_intIOR;
            n = Vector3f(0, 0, 1);
        }
        float costhe1, costhe2;
        Vector3f ome_t;
        float rate_c = pow(eta1 / eta2, 2.0f);
        rate_c *= 1.0f - pow(bRec.wi.dot(n), 2.0f);
        rate_c = 1.0f - rate_c;
        bRec.wo = Vector3f(
            -bRec.wi.x(),
            -bRec.wi.y(),
             bRec.wi.z()
        );
        bRec.eta = 1.0f;
        if (rate_c > 0.0f) {
            ome_t = - n * sqrt(rate_c);
            ome_t -= (eta1 / eta2) * (bRec.wi - bRec.wi.dot(n) * n);
            costhe1 = fabs(Frame::cosTheta(bRec.wi));
            costhe2 = fabs(Frame::cosTheta(ome_t));
            float rho1, rho2, Fr;
            rho1 = (eta2 * costhe1 - eta1 * costhe2) / (eta2 * costhe1 + eta1 * costhe2);
            rho2 = (eta1 * costhe1 - eta2 * costhe2) / (eta1 * costhe1 + eta2 * costhe2);
            Fr = 0.5f * (rho1 * rho1 + rho2 * rho2);
            if (sample.x() > Fr) {
                bRec.wo = ome_t;
                bRec.eta = eta1 / eta2;
            }
        }
        return Color3f(1.0f);
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
