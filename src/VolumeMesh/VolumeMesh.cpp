//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "VolumeMesh.h"
#include "VolumeMeshIO.h"

//=============================================================================

VolumeMesh::VolumeMesh() = default;

//-----------------------------------------------------------------------------

VolumeMesh::~VolumeMesh() = default;

//-----------------------------------------------------------------------------

bool VolumeMesh::read(const std::string &filename)
{
    VolumeMeshIO reader(filename);
    return reader.read(*this);
}

//-----------------------------------------------------------------------------

void VolumeMesh::update_bounding_sphere()
{
    Vec3f v;
    Vec3f min(std::numeric_limits<float>::max());
    Vec3f max(std::numeric_limits<float>::min());

    for (auto bv_it = bv_iter(); bv_it.valid(); ++bv_it)
    {
        v = vertex(*bv_it);
        min.minimize(v);
        max.maximize(v);
    }

    Vec3f center, diameter;
    float radius;

    center = 0.5f * (min + max);
    diameter = max - min;
    radius = 0.5f * diameter.norm();

    bounding_sphere_ = BoundingSphere(center, radius);
}
