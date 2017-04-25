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
#include <nori/object.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class RoughDielectric : public BSDF {
public:
    RoughDielectric(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

    }

    float DW(const Vector3f &w) const{
        if (Frame::cosTheta(w) <= 0.0f)
            return 0.0f;
        float ct, ct2, st2, al2;
        ct = Frame::cosTheta(w);
        ct2 = ct * ct;
        st2 = 1 - ct2;
        al2 = m_alpha * m_alpha;
        return exp(-st2 / ct2 / al2) * INV_PI / al2 / ct2 / ct2;

    }

    float G1(const Vector3f &w, const Vector3f &wh) const {
        if (w.dot(wh) * Frame::cosTheta(w) <= 0.0f)
            return 0.0f;

        float b, ct, st;
        ct = Frame::cosTheta(w);
        st = Frame::sinTheta(w);
        if (st <= 0.0f)
            b = 100.0f;
        else
            b = ct / m_alpha / st;
        float b2 = b * b;
        if (b < 1.6f)
            return (3.535f * b + 2.181f * b2) / (1.0f + 2.276f * b + 2.577f * b2);
        else
            return 1.0f;
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {

        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0.0f)
            return Color3f(0.0f);

        Color3f res;

        Vector3f wo = bRec.wo;
        if (wo.z() < 0.f)
            wo.z() = -wo.z();

        Vector3f wh = bRec.wi + wo;
        wh.normalize();

        res = Color3f(1.0f) * DW(wh);
        res *= fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
        res *= G1(bRec.wi, wh) * G1(wo, wh);
        res /= 4.0f * Frame::cosTheta(bRec.wi) * Frame::cosTheta(wo) * Frame::cosTheta(wh);
        res /= 2.0f;
        // todo: Frame::cosTheta(wh)

        return res;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0.0f)
            return 0.0f;

        Vector3f wo = bRec.wo;
        if (wo.z() < 0.f)
            wo.z() = -wo.z();

        float res = 1.0f;
        Vector3f wh = bRec.wi + wo;
        wh.normalize();
        res *= DW(wh);
        res *= Frame::cosTheta(wh);
        res /= 4.0f * wh.dot(wo);
        res /= 2.0f;
        return res;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        if (Frame::cosTheta(bRec.wi) <= 0.0f)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;

        Point2f new_sample = Point2f(_sample.x(), _sample.y());

        Vector3f wh = Warp::squareToBeckmann(new_sample, m_alpha);
        if (Warp::squareToBeckmannPdf(wh, m_alpha) <= 0.0f)
            return Color3f(0.0f);
        if (G1(bRec.wi, wh) <= 0.0f)
            return Color3f(0.0f);
        float sample2 = bRec.sampler->next1D();

        bRec.wo = 2 * wh * bRec.wi.dot(wh) - bRec.wi;
        if (bRec.wo.z() <= 0.0f)
            return Color3f(0.0f);

        if (sample2 < 0.5f)
            bRec.wo.z() = -bRec.wo.z();
        return Color3f(1.0f);


        float eta1, eta2;
        Vector3f n;
        if (Frame::cosTheta(bRec.wi) <= 0.0f) {
            eta1 = m_intIOR;
            eta2 = m_extIOR;
        }
        else {
            eta1 = m_extIOR;
            eta2 = m_intIOR;
        }

        float F = fresnel(fabs(wh.dot(bRec.wi)), eta1, eta2);

        if (sample2 > F) {
            float c = bRec.wi.dot(wh);
            float eta = eta1 / eta2;
            float weight0 = sqrt(1 + eta * (c * c - 1));
            if (Frame::cosTheta(bRec.wi) <= 0.0f)
                weight0 = -weight0;
            weight0 = eta * c - weight0;
            bRec.wo = weight0 * wh - eta * bRec.wi;
            bRec.eta = eta1 / eta2;
        }
        else {
            bRec.wo = 2 * wh * bRec.wi.dot(wh) - bRec.wi;
            bRec.eta = 1.0f;
        }

        if (G1(bRec.wo, wh) <= 0.0f)
            return Color3f(0.0f);

        return Color3f(1.0f);

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        double pdf_rec = pdf(bRec);
        if (pdf_rec <= 0.0f)
            return Color3f(0.0f);
        else
            return eval(bRec) * fabs(Frame::cosTheta(bRec.wo)) / pdf_rec;
    }

    bool isDiffuse() const {
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "RoughDielectric[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(RoughDielectric, "roughdielectric");
NORI_NAMESPACE_END
