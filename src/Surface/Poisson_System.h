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

double solve_poisson_system(SurfaceMesh &mesh, int minpoint, int function);

double poisson_function(Point &p, int function);

double laplace_of_poisson_function(Point &p, int function);

double franke_function(double x, double y);

double laplace_franke_function(double x, double y);

double spherical_harmonic_function(double x, double y, double z);

double spherical_harmonic_function_scaled(double x, double y, double z);

double inverse_mean_edgelenth(SurfaceMesh &mesh);

void solve_laplace_equation(SurfaceMesh &mesh, int face_point);
