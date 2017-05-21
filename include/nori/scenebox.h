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
#include <nori/emitter.h>
#include <nori/intersection.h>

NORI_NAMESPACE_BEGIN

struct Intersection;

class SceneBox : public NoriObject {
public:

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return ESceneBox; }

    /// Return a pointer to an attached area emitter instance (const version)
    virtual const Emitter *getEmitter() const = 0;

    /// Return a pointer to the BSDF associated with this mesh
    virtual const BSDF *getBSDF() const = 0;

    virtual bool rayIntersect(const Ray3f &ray, Intersection &its) const = 0;

    virtual void activate() = 0;

};

NORI_NAMESPACE_END
