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

#pragma once

#include <map>
#include <nori/object.h>
#include <nori/sampler.h>
#include <nori/volumemedia.h>

NORI_NAMESPACE_BEGIN

class VolumeSurface : public NoriObject {
public:
    VolumeSurface(const PropertyList &props);
    EClassType getClassType() const { return EVolumeSurface; }
    std::string toString() const;
    void activate(std::map<std::string, const VolumeMedia*> &volumeMedia);
    const VolumeMedia *dirMedia(Vector3f wi) const;

private:
    std::string m_intMediaName, m_extMediaName;
    const VolumeMedia *m_intMedia;
    const VolumeMedia *m_extMedia;
};

NORI_NAMESPACE_END
