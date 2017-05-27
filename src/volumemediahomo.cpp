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

#include <nori/volumemedia.h>

NORI_NAMESPACE_BEGIN

class VolumeMediaHomo : public VolumeMedia {
public:
    VolumeMediaHomo(const PropertyList &props) {
        m_name = props.getString("name", "");
        m_bsdf = (BSDF *)(NoriObjectFactory::createInstance("volumeHG", props));
        theta_t = props.getFloat("theta_t", 0.1f);
    }

    ~VolumeMediaHomo() {
        delete m_bsdf;
    }

    std::string getName() const{
        return m_name;
    }

    float sample_dis(Sampler* sampler) const {
        if (theta_t > Epsilon) {
            float dis = sampler->next1D();
            dis = -(log(dis) / theta_t);
            return dis;
        }
        else
            return std::numeric_limits<float>::infinity();
    }

    BSDF* getBSDF(Point3f &p) const {
        return m_bsdf;
    }

    std::string toString() const {
        return tfm::format(
            "VolumeMedia Homogeneous[\n"
            " theta_t = %f\n"
            "]",
            theta_t
        );
    }

private:
    std::string m_name;
    BSDF* m_bsdf = nullptr;
    float theta_t;
};

NORI_REGISTER_CLASS(VolumeMediaHomo, "volumemediahomo");
NORI_NAMESPACE_END
