//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================
#include <VolumeMesh.h>
//=============================================================================

// https://www.dropbox.com/s/dzz7je2cbclq5gy/LB3D.pdf?dl=1 section 8.2
double solve_eigenvalue_problem(VolumeMesh &mesh_, Eigen::VectorXd &evalues,
                                int face_point, int cell_point,
                                std::string meshname = "default");

void analytic_eigenvalues_unitBall(Eigen::VectorXd &eval, int n);
