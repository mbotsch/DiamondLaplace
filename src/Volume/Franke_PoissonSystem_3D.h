//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include "VolumeMesh.h"

//=============================================================================

double franke3d(VolumeMesh::PointT vec);

double laplace_franke3d(VolumeMesh::PointT vec);

double solve_franke_poisson(VolumeMesh& mesh_, Eigen::VectorXd& result,
                            int face_point, int cell_point);

void solve_laplace_equation(VolumeMesh& mesh_, int face_point, int cell_point);