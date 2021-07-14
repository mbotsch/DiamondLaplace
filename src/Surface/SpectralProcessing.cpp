//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#include "SpectralProcessing.h"
#include "LaplaceConstruction.h"
#include <iostream>
#include <Eigen/Sparse>
//=============================================================================

using namespace pmp;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

enum InsertedPoint
{
    Centroid = 0,
    AreaMinimizer = 2
};

void SpectralProcessing::compute_eigen_vectors(unsigned int min_point)
{
    if (!mesh.n_vertices())
        return;
    setup_stiffness_matrix(mesh, L_, min_point);
    setup_mass_matrix(mesh, M_, min_point);
}

//----------------------------------------------------------------------------

double SpectralProcessing::analyze_eigenvectors()
{
    auto points = mesh.vertex_property<Point>("v:point");

    // comparing eigenvectors up to the 9th Band of legendre polynomials

    int band = 8;

    double error;
    double sum = 0.0;

    Eigen::VectorXd y(mesh.n_vertices());
    Eigen::SparseMatrix<double> L, M;
    M = M_;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);

    for (int l = 1; l <= band; l++)
    {
        double eval = -l * (l + 1);
        for (int m = -l; m <= l; m++)
        {
            for (auto v : mesh.vertices())
            {
                y(v.idx()) = sphericalHarmonic(points[v], l, m);
            }
            y.normalize();
            Eigen::MatrixXd X = solver.solve(L_ * y);
            error =
                (y - 1.0 / eval * X).transpose() * M_ * (y - 1.0 / eval * X);
            sum += error;
        }
    }

    if (min_point_ == AreaMinimizer)
    {
        std::cout << "Error SH (our point): " << sum << std::endl;
    }
    else
    {
        std::cout << "Error SH (centroid): " << sum << std::endl;
    }

    return sum;
}
//----------------------------------------------------------------------------

double factorial(int n)
{
    if (n == 0)
        return 1.0;
    return (double)n * factorial(n - 1);
}

//----------------------------------------------------------------------------

double SpectralProcessing::scale(int l, int m)
{
    double temp =
        ((2.0 * l + 1.0) * factorial(l - m)) / (4.0 * M_PI * factorial(l + m));
    return sqrt(temp);
}

//----------------------------------------------------------------------------

double SpectralProcessing::P(int l, int m, double x)
{
    // evaluate an Associated Legendre Polynomial P(l,m,x) at x
    double pmm = 1.0;
    if (m > 0)
    {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= m; i++)
        {
            pmm *= (-fact) * somx2;
            fact += 2.0;
        }
    }
    if (l == m)
        return pmm;
    double pmmp1 = x * (2.0 * m + 1.0) * pmm;
    if (l == m + 1)
        return pmmp1;
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ++ll)
    {
        pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

//----------------------------------------------------------------------------

double SpectralProcessing::sphericalHarmonic(Point p, int l, int m)
{

    // l is the band, range [0..n]
    // m in the range [-l..l]
    // transform cartesian to spherical coordinates, assuming r = 1

    double phi = atan2(p[0], p[2]) + M_PI;
    double theta = acos(p[1] / norm(p));

    // return a point sample of a Spherical Harmonic basis function
    const double sqrt2 = sqrt(2.0);
    if (m == 0)
        return scale(l, 0) * P(l, m, cos(theta));
    else if (m > 0)
        return sqrt2 * scale(l, m) * cos(m * phi) * P(l, m, cos(theta));
    else
        return sqrt2 * scale(l, -m) * sin(-m * phi) * P(l, -m, cos(theta));
}

//=============================================================================
