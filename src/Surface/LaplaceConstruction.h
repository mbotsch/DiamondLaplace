//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>

//=============================================================================

using namespace pmp;

//=============================================================================

void setup_stiffness_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S,
                            int minpoint = 0);

void setup_mass_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M,
                       int minpoint = 0);
