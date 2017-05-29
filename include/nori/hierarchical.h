#pragma once

#include <nori/object.h>
#include <nori/common.h>
#include <nori/bitmap.h>

NORI_NAMESPACE_BEGIN

class HierarchicalSampler {
public:
    HierarchicalSampler ();
    HierarchicalSampler (float max_light);
    void setImage(const std::string &file_name);
    float squareToHierarchicalPdf(const Point2f &p);
    Point2f squareToHierarchical(Point2f sample);
    Color3f hit(Point2f sample);
    float get_mean();
private:
    Point2f squareToHierarchical(Point2f sample, int x, int y, int layer);
    typedef Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> hieImage;
    float mean_light;
    std::vector<hieImage> layers;
    hieImage for_pdf;
    Bitmap bitmap;
    float m_max;
};

NORI_NAMESPACE_END
