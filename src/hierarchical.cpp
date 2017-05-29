#include <nori/hierarchical.h>

NORI_NAMESPACE_BEGIN

HierarchicalSampler::HierarchicalSampler () {
    m_max = 100.0f;
    layers.clear();
}

HierarchicalSampler::HierarchicalSampler (float max_light) {
    m_max = max_light;
    layers.clear();
}

void HierarchicalSampler::setImage(const std::string &file_name) {
    bitmap = Bitmap(file_name);
    layers.clear();
    layers.push_back(hieImage());
    layers[0].resize(bitmap.rows(), bitmap.cols());
    for (int y = 0; y < bitmap.rows(); y++)
        for (int x = 0; x < bitmap.cols(); x++){
            float color = 0;
            Color3f c = bitmap.coeff(y, x);
            color += c.x();
            color += c.y();
            color += c.z();
            if (color > m_max) {
                bitmap.coeffRef(y, x) *= m_max / color;
                color = m_max;
            }
            layers[0].coeffRef(y,x) = color;
        }
    mean_light = layers[0].mean();
    for (int y = 0; y < bitmap.rows(); y++)
        for (int x = 0; x < bitmap.cols(); x++)
            layers[0].coeffRef(y, x) /= mean_light;
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
                layers[new_index].coeffRef(y, x) = layers[index].block(y * 2, x * 2, 2, 2).mean();
            }
        index++;
    }
}

Color3f HierarchicalSampler::hit(Point2f sample) {
    float sx = sample.x();
    float sy = 1.0f - sample.y();
    if ((sx < 0.0f) || (sx > 1.0f) || (sy < 0.0f) || (sy > 1.0f))
        return Color3f(0.0f);
    int x, y;
    x = sx * bitmap.cols();
    x = (x > bitmap.cols() - 1)? bitmap.cols() - 1 : x;
    x = (x < 0)? 0 : x;
    y = sy * bitmap.rows();
    y = (y > bitmap.rows() - 1)? bitmap.rows() - 1 : y;
    y = (y < 0)? 0 : y;
    return bitmap.coeff(y, x);
}

float HierarchicalSampler::get_mean() {
    return mean_light;
}

float HierarchicalSampler::squareToHierarchicalPdf(const Point2f &p) {
    if (layers[0].rows() == 0 || layers[0].cols() == 0)
        return 1;
    int row, col;
    row = (1 - p.y()) * layers[0].rows();
    col = p.x() * layers[0].cols();
    if (row >= layers[0].rows())
        row = layers[0].rows() - 1;
    if (col >= layers[0].cols())
        col = layers[0].cols() - 1;
    return layers[0].coeff(row, col);
}

Point2f HierarchicalSampler::squareToHierarchical(Point2f sample, int x, int y, int layer) {
    float sx = sample.x();
    float sy = sample.y();
    sx = (sx < 0) ? 0 : sx;
    sx = (sx > 1.0f) ? 1.0f : sx;
    sy = (sy < 0) ? 0 : sy;
    sy = (sy > 1.0f) ? 1.0f : sy;
    if (layer < 0)
        return Point2f(sx, sy);
    float c00, c01, c10, c11;
    c00 = layers[layer].coeff(y, x);
    c01 = layers[layer].coeff(y, x + 1);
    c10 = layers[layer].coeff(y + 1, x);
    c11 = layers[layer].coeff(y + 1, x + 1);
    float t0, t1;
    t0 = (c00 + c10) / (c00 + c01 + c10 + c11);
    float start_x, start_y;
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

Point2f HierarchicalSampler::squareToHierarchical(Point2f sample) {
    Point2f res = squareToHierarchical(sample, 0, 0, layers.size() - 1);
    return (Point2f(res.x(), 1.0f - res.y()));
}

NORI_NAMESPACE_END