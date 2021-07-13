//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>

//=============================================================================

using namespace pmp;

//=============================================================================

//!compute the transformation matrix for an arbitrary mesh that inserts a chosen point per face .
void setup_prolongation_matrix(SurfaceMesh &mesh,
                               Eigen::SparseMatrix<double> &P);

//----------------------------------area computations------------------------------------------------------

//! barycenter/centroid of mesh, computed as area-weighted mean of vertices.
Point my_centroid(const SurfaceMesh &mesh);

//! computes the area of a face.
double face_area(const SurfaceMesh &mesh, Face f);

//! surface area of the mesh.
Scalar my_surface_area(const SurfaceMesh &mesh);

//! Computes the squared triangle area minimizing points and its convex combination weights
//! for each face and stores it in a prior defined property.
void setup_face_point_properties(SurfaceMesh &mesh, unsigned int min_point);

//------------------------Point and Weight minimizer -----------------------------------------------------------------

void find_area_minimizer_weights(const Eigen::MatrixXd &poly,
                                 Eigen::VectorXd &weights);

//=============================================================================