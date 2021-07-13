//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include <Eigen/Sparse>
#include "../VolumeMesh/VolumeMesh.h"

//=============================================================================

class GeodesicsInHeat_3d
{
public:
    GeodesicsInHeat_3d(VolumeMesh &mesh, unsigned int face_point,
                       unsigned int cell_point);

    ~GeodesicsInHeat_3d();

    double compute_geodesics(const int vertex, Eigen::VectorXd &dist,
                             Eigen::VectorXd &orthodist, bool max_e_len = true);

    void distance_to_texture_coordinates() const;
private:
    VolumeMesh &mesh_;

    unsigned int face_point_, cell_point_;

    double maxEdgeLength(const VolumeMesh &mesh);

    double meanEdgeLength(const VolumeMesh &mesh);
};
