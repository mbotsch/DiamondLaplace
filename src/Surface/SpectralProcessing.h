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

class SpectralProcessing
{
public:
    SpectralProcessing(SurfaceMesh& mesh)
        : mesh(mesh), vertices_(0), faces_(0), min_point_(0)
    {
    }

    void compute_eigen_vectors(unsigned int min_point);

    void update_eigenvectors(unsigned int min_point)
    {
        if (mesh.n_faces() != faces_ || mesh.n_vertices() != vertices_ ||
            min_point != min_point_)
        {
            min_point_ = min_point;
            vertices_ = mesh.n_vertices();
            faces_ = mesh.n_faces();
            compute_eigen_vectors(min_point);
        }
    }

    //// renormalisation constant for SH function
    double scale(int l, int m);

    // spherical harmonic function at a point p for a band l and its range m
    // see: http://silviojemma.com/public/papers/lighting/spherical-harmonic-lighting.pdf
    double sphericalHarmonic(Point p, int l, int m);

    // evaluate an Associated Legendre Polynomial P(l,m,x) at x
    double P(int l, int m, double x);

    double analyze_eigenvectors();

private:
    SurfaceMesh& mesh;

    Eigen::SparseMatrix<double> M_, L_;

    unsigned int vertices_, faces_, min_point_;
};

//=============================================================================
