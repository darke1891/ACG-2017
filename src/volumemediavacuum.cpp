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

class VolumeMediaVacuum : public VolumeMedia {
public:
    VolumeMediaVacuum(const PropertyList &props) {
        m_name = "";
    }

    std::string getName() const{
        return m_name;
    }

    float sample_dis(Sampler* sampler) const {
        return std::numeric_limits<float>::infinity();
    }

    BSDF* getBSDF(Point3f &p) const {
        return nullptr;
    }

    std::string toString() const {
        return tfm::format(
            "VolumeMedia Vacuum[\n"
            "]"
        );
    }

private:
    std::string m_name;
};

NORI_REGISTER_CLASS(VolumeMediaVacuum, "volumemediavacuum");
NORI_NAMESPACE_END
