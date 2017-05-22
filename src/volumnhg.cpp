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

class VolumnHG : public BSDF {
public:
    VolumnHG(const PropertyList &propList) {
        m_g = propList.getFloat("g", 0.0f);
        float almost_one = 0.9999f;
        m_g = (m_g > almost_one)? almost_one : m_g;
        m_g = (m_g < -almost_one)? -almost_one : m_g;
    }

    float HG_pdf(const BSDFQueryRecord &bRec) const {
        float costheta = -Frame::cosTheta(bRec.wo);
        float res = 1.0f - m_g * m_g;
        res /= 4.0f * pow(1.0f + m_g * m_g - 2 * m_g * costheta, 1.5f);
        res *= INV_PI;
        return res;
    }

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) != 1.0f)
            return Color3f(0.0f);

        /* The BRDF is simply the albedo / pi */
        return Color3f(1.0f) * HG_pdf(bRec);
    }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        /* This is a smooth BRDF -- return zero if the measure
           is wrong, or when queried for illumination on the backside */
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) != 1.0f)
            return 0.0f;


        /* Importance sampling density wrt. solid angles:
           cos(theta) / pi.

           Note that the directions in 'bRec' are in local coordinates,
           so Frame::cosTheta() actually just returns the 'z' component.
        */
        return HG_pdf(bRec);
    }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        if (Frame::cosTheta(bRec.wi) != 1.0f)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        /* Warp a uniformly distributed sample on [0,1]^2
           to a direction on a cosine-weighted hemisphere */
        if (m_g == 0.0f) {
            bRec.wo = Warp::squareToUniformSphere(sample);
        }
        else {
            float tmp = 1.0f - m_g * m_g;
            float h, r;
            float theta = 2.0f * M_PI * sample.x();

            tmp /= 1.0f - m_g + 2.0f * m_g * sample.y();
            tmp = 1.0f + m_g * m_g - tmp * tmp;
            h = tmp / 2.0f / m_g;
            r = sqrt(1.0f - h * h);
            bRec.wo = Vector3f(r * cos(theta), r * sin(theta), -h);
        }


        /* Relative index of refraction: no change */
        bRec.eta = 1.0f;

        /* eval() / pdf() * cos(theta) = albedo. There
           is no need to call these functions. */
        float pdf_rec = pdf(bRec);
        if (pdf_rec <= 0.0f)
            return Color3f(0.0f);
        else
            return eval(bRec) / pdf_rec;
    }

    bool isDiffuse() const {
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "VolumnHG[\n"
            " g = %f\n"
            "]",
            m_g
        );
    }
private:
    float m_g;
};

NORI_REGISTER_CLASS(VolumnHG, "volumnHG");
NORI_NAMESPACE_END
