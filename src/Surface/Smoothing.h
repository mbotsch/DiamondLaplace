//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "LaplaceConstruction.h"

//=============================================================================

using namespace pmp;

//=============================================================================

class Smoothing
{
public:
    Smoothing(SurfaceMesh& mesh)
        : mesh_(mesh), vertices_(0), faces_(0), min_point_(0)
    {
    }

    void implicit_smoothing(Scalar timestep, unsigned int min_point);

private:
    void update_Laplace(unsigned int min_point)
    {
        if (mesh_.n_faces() != faces_ || mesh_.n_vertices() != vertices_ ||
            min_point_ != min_point)
        {
            vertices_ = mesh_.n_vertices();
            faces_ = mesh_.n_faces();
            min_point_ = min_point;

            std::cout << "Stiffness matrix has been updated" << std::endl;
            setup_stiffness_matrix(mesh_, S_, min_point_);
        }
    }

private:
    SurfaceMesh& mesh_;

    Eigen::SparseMatrix<double> S_;
    unsigned int vertices_, faces_, min_point_;
};

//=============================================================================
