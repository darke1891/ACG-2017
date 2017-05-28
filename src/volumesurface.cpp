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

#include <nori/volumesurface.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

VolumeSurface::VolumeSurface(const PropertyList &props) {
    m_intMediaName = props.getString("intMedia", "");
    m_extMediaName = props.getString("extMedia", "");
}

void VolumeSurface::activate(std::map<std::string, const VolumeMedia*> &volumeMedia) {
    auto its = volumeMedia.find(m_intMediaName);
    if (its == volumeMedia.end())
        throw NoriException("No integrator was specified!");
    m_intMedia = its->second;
    its = volumeMedia.find(m_extMediaName);
    if (its == volumeMedia.end())
        throw NoriException("No integrator was specified!");
    m_extMedia = its->second;
}

const VolumeMedia *VolumeSurface::dirMedia(Vector3f wi) const {
    if (Frame::cosTheta(wi) <= 0.0f)
        return m_extMedia;
    else
        return m_intMedia;
}

std::string VolumeSurface::toString() const {
    return tfm::format(
        "Volume Surface[\n"
            "  intMedia = \"%s\",\n"
            "  extMedia = %s,\n"
        "]",
        m_intMediaName,
        m_extMediaName
    );
}

NORI_REGISTER_CLASS(VolumeSurface, "volumesurfaceN");
NORI_NAMESPACE_END
