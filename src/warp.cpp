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

void Warp::HierarchicalSampler::setImage(Bitmap &bitmap) {
    layers.clear();
    layers.push_back(hieImage());
    layers[0].resize(bitmap.rows(), bitmap.cols());
    float sum_color = 0;
    for (int y = 0; y < bitmap.rows(); y++)
        for (int x = 0; x < bitmap.cols(); x++){
            float color = 0;
            Color3f c = bitmap.coeff(y, x);
            color += c.x();
            color += c.y();
            color += c.z();
            layers[0].coeffRef(y,x) = color;
            sum_color += layers[0].coeff(y,x);
        }
    sum_color /= bitmap.rows() * bitmap.cols();
    for (int y = 0; y < bitmap.rows(); y++)
        for (int x = 0; x < bitmap.cols(); x++)
            layers[0].coeffRef(y, x) /= sum_color;
    int index = 0;
    while (true) {
        int new_index = index + 1;
        int rows = layers[index].rows();
        int cols = layers[index].cols();
        if (rows <= 2)
            break;
        layers.push_back(hieImage());
        int new_rows = rows / 2;
        int new_cols = cols / 2;
        layers[new_index].resize(new_rows, new_cols);
        for (int y=0; y<new_rows; y++)
            for (int x=0; x<new_cols; x++) {
                double new_color = 0;
                new_color += layers[index].coeff(y * 2, x * 2);
                new_color += layers[index].coeff(y * 2 + 1, x * 2);
                new_color += layers[index].coeff(y * 2, x * 2 + 1);
                new_color += layers[index].coeff(y * 2 + 1, x * 2 + 1);
                new_color /= 4;
                layers[new_index].coeffRef(y, x) = new_color;
            }
        index++;
    }
}

float Warp::HierarchicalSampler::squareToHierarchicalPdf(const Point2f &p) {
    if (for_pdf.rows() == 0 || for_pdf.cols() == 0)
        return 1;
    int row, col;
    row = (1 - p.y()) * for_pdf.rows();
    col = p.x() * for_pdf.cols();
    if (row >= for_pdf.rows())
        row = for_pdf.rows() - 1;
    if (col >= for_pdf.cols())
        col = for_pdf.cols() - 1;
    return for_pdf.coeff(row, col);
}

void Warp::HierarchicalSampler::setTestLayer(int xres, int yres) {
    if (layers.size() == 0)
        return;
    
    int index = 0;
    while (layers[index].rows() > yres)
        index++;
    for_pdf = layers[index];
}

Point2f Warp::HierarchicalSampler::squareToHierarchical(Point2f sample, int x, int y, int layer) {
    double sx = sample.x();
    double sy = sample.y();
    sx = (sx < 0) ? 0 : sx;
    sx = (sx > 1.0f) ? 1.0f : sx;
    sy = (sy < 0) ? 0 : sy;
    sy = (sy > 1.0f) ? 1.0f : sy;
    if (layer < 0)
        return Point2f(sx, sy);
    double c00, c01, c10, c11;
    c00 = layers[layer].coeff(y, x);
    c01 = layers[layer].coeff(y, x + 1);
    c10 = layers[layer].coeff(y + 1, x);
    c11 = layers[layer].coeff(y + 1, x + 1);
    double t0, t1;
    t0 = (c00 + c10) / (c00 + c01 + c10 + c11);
    double start_x, start_y;
    if (sx < t0) {
        sx /= t0;
        t1 = c00 / (c00 + c10);
        x = x * 2;
        start_x = 0;
    }
    else {
        sx = (sx - t0) / (1 - t0);
        t1 = c01 / (c01 + c11);
        x = x * 2 + 2;
        start_x = 0.5f;
    }
    if (sy < t1) {
        sy /= t1;
        y = y * 2;
        start_y = 0;
    }
    else {
        sy = (sy - t1) / (1 - t1);
        y = y * 2 + 2;
        start_y = 0.5f;
    }
    Point2f next_p = squareToHierarchical(Point2f(sx, sy), x, y, layer-1);
    Point2f res = (Point2f(next_p.x() * 0.5f + start_x, next_p.y() * 0.5f + start_y));
    return res;

}

Point2f Warp::HierarchicalSampler::squareToHierarchical(Point2f sample) {
    Point2f res = squareToHierarchical(sample, 0, 0, layers.size() - 1);
    return (Point2f(res.x(), 1.0f - res.y()));
}

Point2f Warp::squareToHierarchical(const Point2f &sample, HierarchicalSampler &h_sampler) {
    Point2f s = Point2f(sample.x(), sample.y());
    return h_sampler.squareToHierarchical(s);
}

float Warp::squareToHierarchicalPdf(const Point2f &p, HierarchicalSampler &h_sampler) {
    Point2f p2 = Point2f(p.x(), p.y());
    return h_sampler.squareToHierarchicalPdf(p2);
}

Warp::HierarchicalSampler::HierarchicalSampler () {
    layers.clear();
}

NORI_NAMESPACE_END
