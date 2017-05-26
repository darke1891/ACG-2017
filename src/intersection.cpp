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

#include <nori/intersection.h>

NORI_NAMESPACE_BEGIN

/// Return a pointer to an attached area emitter instance (const version)
const Emitter *Intersection::getEmitter() const {
    if (mesh)
        return mesh->getEmitter();
    else if (sbox)
        return sbox->getEmitter();
    else
        return m_emitter;
}

/// Return a pointer to the BSDF associated with this mesh
const BSDF *Intersection::getBSDF() const {
    if (mesh)
        return mesh->getBSDF();
    else if (sbox)
        return sbox->getBSDF();
    else
        return m_bsdf;
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

NORI_NAMESPACE_END
