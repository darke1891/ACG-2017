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

#include <nori/object.h>
#include <nori/sampler.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/ray.h>

NORI_NAMESPACE_BEGIN

class VolumeMedia : public NoriObject {
public:
    EClassType getClassType() const { return EVolumeMedia; }
    virtual std::string getName() const = 0;
    virtual float sample_dis(const Ray3f &ray, float t, Sampler* sampler) const = 0;
    virtual BSDF* getBSDF() const = 0;
    virtual Emitter* getEmitter() const = 0;
};

NORI_NAMESPACE_END
