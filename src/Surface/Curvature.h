//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"

//=============================================================================

using namespace pmp;

//=============================================================================

class Curvature
{
public:
    Curvature(SurfaceMesh& mesh, bool compare)
        : mesh_(mesh), compare_to_sphere(compare)
    {
    }

    //! Visualizes the mean curvature of our mesh.
    void visualize_curvature(unsigned int min_point_);

    double compute_curvature_error(unsigned int min_point_);

private:
    SurfaceMesh& mesh_;
    bool compare_to_sphere;

    //! convert curvature values ("v:curv") to 1D texture coordinates
    void curvature_to_texture_coordinates() const;
};

//=============================================================================
