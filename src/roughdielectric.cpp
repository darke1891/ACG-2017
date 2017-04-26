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

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    float DW(const Vector3f &w, float alpha) const{
        if (Frame::cosTheta(w) <= 0.0f)
            return 0.0f;
        float ct, ct2, st2, al2;
        ct = Frame::cosTheta(w);
        ct2 = ct * ct;
        st2 = 1 - ct2;
        al2 = alpha * alpha;
        return exp(-st2 / ct2 / al2) * INV_PI / al2 / ct2 / ct2;

    }

    float G1(const Vector3f &w, const Vector3f &wh, float alpha) const {
        if (w.dot(wh) * Frame::cosTheta(w) <= 0.0f)
            return 0.0f;

        float b, ct, st;
        ct = fabs(Frame::cosTheta(w));
        st = Frame::sinTheta(w);
        if (st <= 0.0f)
            b = 100.0f;
        else
            b = ct / alpha / st;
        float b2 = b * b;
        if (b < 1.6f)
            return (3.535f * b + 2.181f * b2) / (1.0f + 2.276f * b + 2.577f * b2);
        else
            return 1.0f;
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {

        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) == 0.0f
            || Frame::cosTheta(bRec.wo) == 0.0f)
            return Color3f(0.0f);

        float alpha = (1.2f - 0.2f * sqrt(fabs(Frame::cosTheta(bRec.wi)))) * m_alpha;

        Vector3f wi, wo, n, wh;
        float eta1, eta2;
        wi = bRec.wi;
        wo = bRec.wo;
        if (Frame::cosTheta(wi) > 0.0f) {
            eta1 = m_extIOR;
            eta2 = m_intIOR;
            n = Vector3f(0, 0, 1.0f);
        }
        else {
            eta1 = m_intIOR;
            eta2 = m_extIOR;
            n = Vector3f(0, 0, -1.0f);
        }
        float res;
        if (Frame::cosTheta(wi) * Frame::cosTheta(wo) > 0.0f) {
            wh = wi + wo;
            wh.normalize();
            if (Frame::cosTheta(wi) < 0.0f)
                wh *= -1.0f;
            if (G1(wi, wh, alpha) <= 0.0f)
                return Color3f(0.0f);
            if (G1(wo, wh, alpha) <= 0.0f)
                return Color3f(0.0f);
            if (Frame::cosTheta(wh) <= 0.0f)
                return Color3f(0.0f);

            float F = fresnel(fabs(wh.dot(wi)), eta1, eta2);
            res = F;
            res *= DW(wh, alpha);
            res *= G1(wi, wh, alpha) * G1(wo, wh, alpha);
            res /= 4.0f * fabs(Frame::cosTheta(bRec.wi) * Frame::cosTheta(wo) * Frame::cosTheta(wh));

        }
        else {
            wh = -(eta1 * wi + eta2 * wo);
            wh.normalize();

            if (G1(wi, wh, alpha) <= 0.0f)
                return Color3f(0.0f);
            if (G1(wo, wh, alpha) <= 0.0f)
                return Color3f(0.0f);
            if (Frame::cosTheta(wh) <= 0.0f)
                return Color3f(0.0f);

            float F = fresnel(fabs(wh.dot(wi)), eta1, eta2);
            float eta = eta1 / eta2;
            float weight0 = eta * wi.dot(wh) + wo.dot(wh);

            res = 1.0f - F;
            res /= weight0 * weight0;
            res *= DW(wh, alpha);
            res *= G1(wi, wh, alpha) * G1(wo, wh, alpha);
            res /= fabs(Frame::cosTheta(wi) * Frame::cosTheta(wo) * Frame::cosTheta(wh));
            res *= fabs(wo.dot(wh) * wi.dot(wh));
        }

        return Color3f(res);
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) == 0.0f
            || Frame::cosTheta(bRec.wo) == 0.0f)
            return 0.0f;

        float alpha = (1.2f - 0.2f * sqrt(fabs(Frame::cosTheta(bRec.wi)))) * m_alpha;

        Vector3f wi, wo, n, wh;
        float eta1, eta2;
        wi = bRec.wi;
        wo = bRec.wo;
        if (Frame::cosTheta(wi) > 0.0f) {
            eta1 = m_extIOR;
            eta2 = m_intIOR;
            n = Vector3f(0, 0, 1.0f);
        }
        else {
            eta1 = m_intIOR;
            eta2 = m_extIOR;
            n = Vector3f(0, 0, -1.0f);
        }
        float res;
        if (Frame::cosTheta(wi) * Frame::cosTheta(wo) > 0.0f) {
            wh = wi + wo;
            wh.normalize();
            if (Frame::cosTheta(wi) < 0.0f)
                wh *= -1.0f;
            if (G1(wi, wh, alpha) <= 0.0f)
                return 0.0f;
            if (G1(wo, wh, alpha) <= 0.0f)
                return 0.0f;
            if (Frame::cosTheta(wh) <= 0.0f)
                return 0.0f;

            float F = fresnel(fabs(wh.dot(wi)), eta1, eta2);
            res = F;
            res /= 4.0f * fabs(wh.dot(wo));

            res *= DW(wh, alpha);
            res *= fabs(Frame::cosTheta(wh));

        }
        else {
            wh = -(eta1 * wi + eta2 * wo);
            wh.normalize();

            if (G1(wi, wh, alpha) <= 0.0f)
                return 0.0f;
            if (G1(wo, wh, alpha) <= 0.0f)
                return 0.0f;
            if (Frame::cosTheta(wh) <= 0.0f)
                return 0.0f;

            float F = fresnel(fabs(wh.dot(wi)), eta1, eta2);
            res = 1.0f - F;
            float eta = eta1 / eta2;
            float weight0 = eta * wi.dot(wh) + wo.dot(wh);
            res /= weight0 * weight0;
            res *= fabs(wo.dot(wh));

            res *= DW(wh, alpha);
            res *= fabs(Frame::cosTheta(wh));
        }

        return res;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        if (Frame::cosTheta(bRec.wi) == 0.0f)
            return Color3f(0.0f);

        bRec.measure = ESolidAngle;
        float alpha = (1.2f - 0.2f * sqrt(fabs(Frame::cosTheta(bRec.wi)))) * m_alpha;

        Point2f new_sample = Point2f(_sample.x(), _sample.y());

        Vector3f wh = Warp::squareToBeckmann(new_sample, alpha);
        if (Warp::squareToBeckmannPdf(wh, alpha) <= 0.0f)
            return Color3f(0.0f);
        if (G1(bRec.wi, wh, alpha) <= 0.0f)
            return Color3f(0.0f);
        float sample2 = bRec.sampler->next1D();

        float eta1, eta2;
        if (Frame::cosTheta(bRec.wi) < 0.0f) {
            eta1 = m_intIOR;
            eta2 = m_extIOR;
        }
        else {
            eta1 = m_extIOR;
            eta2 = m_intIOR;
        }

        float F = fresnel(fabs(wh.dot(bRec.wi)), eta1, eta2);
        if (sample2 > F) {
            float eta = eta1 / eta2;
            float weight0 = bRec.wi.dot(wh);
            weight0 = eta * eta * (1 - weight0 * weight0);
            weight0 = sqrt(1.0f - weight0);
            if (Frame::cosTheta(bRec.wi) < 0.0f)
                weight0 *= -1.0f;
            Vector3f wt = - weight0 * wh;
            wt -= eta * (bRec.wi - bRec.wi.dot(wh) * wh);
            wt.normalize();

            bRec.wo = wt;
            bRec.eta = eta;
            if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(bRec.wi) >= 0.0f)
                return Color3f(0.0f);
        }
        else {
            bRec.wo = 2 * wh * bRec.wi.dot(wh) - bRec.wi;
            bRec.eta = 1.0f;
            if (Frame::cosTheta(bRec.wo) * Frame::cosTheta(bRec.wi) <= 0.0f)
                return Color3f(0.0f);
        }

        if (G1(bRec.wo, wh, alpha) <= 0.0f)
            return Color3f(0.0f);

        double pdf_rec = pdf(bRec);
        if (pdf_rec <= 0.0f)
            return Color3f(0.0f);

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
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(RoughDielectric, "roughdielectric");
NORI_NAMESPACE_END
