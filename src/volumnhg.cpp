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
        m_color = propList.getColor("color", Color3f(1.0f));
        m_a = propList.getFloat("albedo", 0.9f);
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

    Color3f eval(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) != 1.0f)
            return Color3f(0.0f);

        return Color3f(1.0f) * HG_pdf(bRec) * m_a * m_color;
    }

    float pdf(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) != 1.0f)
            return 0.0f;
        return HG_pdf(bRec);
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        if (Frame::cosTheta(bRec.wi) != 1.0f)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

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


        bRec.eta = 1.0f;

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
    Color3f m_color;
    float m_a;
};

NORI_REGISTER_CLASS(VolumnHG, "volumnHG");
NORI_NAMESPACE_END
