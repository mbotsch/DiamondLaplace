//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include "../VolumeMesh/VolumeMesh.h"

//=============================================================================

void setup_3D_cell_face_prolongation_matrix(VolumeMesh &mesh,
                                            Eigen::SparseMatrix<double> &P,
                                            Eigen::SparseMatrix<double> &Pc,
                                            Eigen::SparseMatrix<double> &Pf);

void setup_3D_stiffness_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &S,
                               int face_point, int cell_point);

void setup_3D_mass_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &M,
                          int face_point, int cell_point);
