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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

float inverseTent(float t) {
    if (t < 0.5f)
        return - 1.0f + sqrt(2 * t);
    else
        return 1.0f - sqrt(2 * (1 - t));
}

Point2f Warp::squareToTent(const Point2f &sample) {
    return Point2f(inverseTent(sample.x()), inverseTent(sample.y()));
}

float Warp::squareToTentPdf(const Point2f &p) {
    if ((p.array() >= -1).all() && (p.array() <= 1).all())
        return (1.0f - fabs(p.x())) * (1.0f - fabs(p.y()));
    else
        return 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float theta = 2 * M_PI * sample.x();
    float r = sqrt(sample.y());
    return Point2f(r * sin(theta), r * cos(theta));
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    return (p.norm() <= 1) ? 1.0f / M_PI : 0.0f;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float theta = 2.0f * M_PI * sample.x();
    float h = sample.y() * 2 - 1.0f;
    float r = sqrt(1.0f - h * h);
    return Vector3f(r * cos(theta), r * sin(theta), h);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return (v.norm() <= 1) ? 0.25f / M_PI : 0.0f;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float theta = 2.0f * M_PI * sample.x();
    float h = sample.y();
    float r = sqrt(1.0f - h * h);
    return Vector3f(r * cos(theta), r * sin(theta), h);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return ((v.z() >= 0) && (v.norm() <= 1)) ? 0.5f / M_PI : 0.0f;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    float h = sample.y();
    h = sqrt(h);
    float theta = 2.0f * M_PI * sample.x();
    float r = sqrt(1.0f - h * h);
    return Vector3f(r * cos(theta), r * sin(theta), h);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    if ((v.z() < 0) || (v.norm() > 1))
        return 0.0f;
    return v.z() / M_PI;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float cth, sth;
    if (sample.y() == 1.0f)
        cth = 0.0f;
    else {
        float t1 = log(1 - sample.y()) * alpha * alpha;
        t1 = 1 - t1;
        cth = sqrt(1.0f / t1);
    }

    float phi = 2.0f * M_PI * sample.x();
    sth = sqrt(1 - cth * cth);
    return Vector3f(sth * cos(phi), sth * sin(phi), cth);
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    if ((m.z() > 0) && (m.norm() <= 1)) {
        float theta = acos(m.z());
        float d = 0.5f / M_PI;
        float t1 = tan(theta) / alpha;
        float t2 = cos(theta);
        d *= 2.0f * exp(-t1 * t1);
        d /= alpha * alpha * t2 * t2 * t2;
        return d;
    }
    else
        return 0.0f;
}

Point2f Warp::squareToHierarchical(const Point2f &sample, HierarchicalSampler &h_sampler) {
    Point2f s = Point2f(sample.x(), sample.y());
    return h_sampler.squareToHierarchical(s);
}

float Warp::squareToHierarchicalPdf(const Point2f &p, HierarchicalSampler &h_sampler) {
    Point2f p2 = Point2f(p.x(), p.y());
    return h_sampler.squareToHierarchicalPdf(p2);
}

NORI_NAMESPACE_END
