//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#include <igl/slice.h>
#include "Poisson_System.h"
//=============================================================================
using namespace std;
//=============================================================================

enum InsertedPoint
{
    Centroid = 0,
    AreaMinimizer = 2
};

enum Function
{

    Franke2d = 0,
    SH = 1
};

//-----------------------------------------------------------------------------

double solve_poisson_system(SurfaceMesh &mesh, int minpoint, int function)
{
    Eigen::SparseMatrix<double> S, M;
    int nv = mesh.n_vertices();
    Eigen::VectorXd b(nv);

    setup_stiffness_matrix(mesh, S, minpoint);
    setup_mass_matrix(mesh, M, minpoint);

    int nb = 0;
    for (auto v : mesh.vertices())
    {
        Point p = mesh.position(v);
        b(v.idx()) = laplace_of_poisson_function(p, function);
        if (mesh.is_boundary(v))
        {
            nb++;
        }
    }

    b = M * b;
    unsigned ins = 0;
    unsigned out = 0;

    Eigen::VectorXi in(mesh.n_vertices() - nb), bound(nb), x(1);

    Eigen::VectorXd solution(mesh.n_vertices());

    x << 0;

    if (nb == 0)
    {
        in.resize(mesh.n_vertices() - 1);
        bound.resize(1);
        for (auto v : mesh.vertices())
        {
            if (v.idx() == mesh.n_vertices() - 1)
            {
                bound(out) = v.idx();
                out++;
                Point p = mesh.position(v);
                b(v.idx()) = poisson_function(p, function);
                solution(v.idx()) = poisson_function(p, function);
            }
            else
            {
                in(ins) = v.idx();
                ins++;
            }
        }
    }
    else
    {
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                in(ins) = v.idx();
                ins++;
            }
            else
            {
                bound(out) = v.idx();
                out++;
                // right-hand side: fix boundary values with franke function of the vertices
                Point p = mesh.position(v);
                b(v.idx()) = poisson_function(p, function);
                solution(v.idx()) = poisson_function(p, function);
            }
        }
    }

    Eigen::SparseMatrix<double> L_in_in, L_in_b;
    Eigen::MatrixXd b_in, b_out;

    igl::slice(S, in, in, L_in_in);
    igl::slice(S, in, bound, L_in_b);
    igl::slice(b, in, x, b_in);
    igl::slice(b, bound, x, b_out);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L_in_in);
    Eigen::MatrixXd X = solver.solve(b_in - L_in_b * b_out);
    int k = 0;
    double error = 0.0;
    if (solver.info() != Eigen::Success)
    {
        std::cout << "Issue: " << solver.info() << std::endl;
        std::cerr << "Could not solve linear system\n";
    }
    else
    {
        // copy solution
        k = 0;
        for (auto v : mesh.vertices())
        {
            if (nb == 0)
            {
                if (v.idx() != mesh.n_vertices() - 1)
                {
                    solution(v.idx()) = X(k);
                    Point p = mesh.position(v);
                    error += pow(X(k) - poisson_function(p, function), 2);
                    k++;
                }
            }
            else
            {
                if (!mesh.is_boundary(v))
                {
                    solution(v.idx()) = X(k);
                    Point p = mesh.position(v);
                    error += pow(X(k) - poisson_function(p, function), 2);
                    k++;
                }
            }
        }
    }
    std::cout << "Error: " << sqrt(error / (double)k) << std::endl;
    return sqrt(error / (double)k);
}

//-----------------------------------------------------------------------------

double poisson_function(Point &p, int function)
{
    switch (function)
    {
        case Franke2d:
            return franke_function(p[0], p[1]);
            break;
        case SH:
            return spherical_harmonic_function_scaled(p[0], p[1], p[2]);
            break;
        default:
            return -1.0;
    }
}

//-----------------------------------------------------------------------------

double laplace_of_poisson_function(Point &p, int function)
{
    switch (function)
    {
        case Franke2d:
            return laplace_franke_function(p[0], p[1]);
            break;
        case SH:
            return spherical_harmonic_function(p[0], p[1], p[2]);
            break;
        default:
            return -1.0;
    }
}

//-----------------------------------------------------------------------------

void sample_weights(SurfaceMesh &mesh, Eigen::MatrixXd &weights, Face face,
                    int n)
{
    int valence = mesh.valence(face);
    Eigen::VectorXd wj(valence);
    weights.resize(n, valence);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < valence; ++j)
        {
            wj(j) = ((double)rand() / (RAND_MAX));
        }
        wj /= wj.sum();
        weights.row(i) = wj;
    }
}

//-----------------------------------------------------------------------------

double franke_function(double x, double y)
{
    double cx2 = (9 * x - 2) * (9 * x - 2);
    double cy2 = (9 * y - 2) * (9 * y - 2);

    double cx1 = (9 * x + 1) * (9 * x + 1);
    double cx7 = (9 * x - 7) * (9 * x - 7);

    double cy3 = (9 * y - 3) * (9 * y - 3);
    double cx4 = (9 * x - 4) * (9 * x - 4);

    double cy7 = (9 * y - 7) * (9 * y - 7);

    return (3. / 4.) * exp(-(1. / 4.) * cx2 - (1. / 4.) * cy2) +
           (3. / 4.) * exp(-(1. / 49.) * cx1 - (9. / 10.) * y - 1. / 10.) +
           (1. / 2.) * exp(-(1. / 4.) * cx7 - (1. / 4.) * cy3) -
           (1. / 5.) * exp(-cx4 - cy7);
}

//-----------------------------------------------------------------------------

double laplace_franke_function(double x, double y)
{

    double wolfram_alpa =
        exp(-(4 - 9 * x) * (4 - 9 * x) - (7 - 9 * y) * (7 - 9 * y)) *
            (-5248.8 * x * x + 4665.6 * x - 5248.8 * y * y + 8164.8 * y -
             4147.2) +
        166.488 * exp(-9 / 4.0 * (9 * x * x - 4 * x + y * (9 * y - 4))) *
            (x * x - 0.444444 * x + y * y - 0.444444 * y + 0.0493827) +
        exp(-1 / 4.0 * (7 - 9 * x) * (7 - 9 * x) -
            9 / 4.0 * (1 - 3 * y) * (1 - 3 * y)) *
            (820.125 * x * x - 1275.75 * x + 820.125 * y * y - 546.75 * y +
             546.75) +
        (7.41771 * x * x + 1.64838 * x - 1.60236) *
            exp(-1 / 49.0 * (9 * x + 1) * (9 * x + 1) - (9 * y) / 10);
    return wolfram_alpa;
}

//-----------------------------------------------------------------------------

double spherical_harmonic_function(double x, double y, double z)
{

    //    y_4,-3 (evalue 20)
    //    return -(3*pow(x, 2.0) - pow(y, 2.0))*y*z; //* (3.0/8.0*sqrt(5.0/M_PI));#

    //y_4,2 (evalue 20)
    //    return -(pow(x, 2.0) - pow(y, 2.0)) * (7.0 * pow(z, 2.0) - 1.0); //* (3.0/8.0*sqrt(5.0/M_PI));#

    // y_3,-1 (evalue 12)
    return -y * (4.0 * pow(z, 2.0) - pow(x, 2.0) - pow(y, 2.0));

    // y_2,0 (evalue 6)
    //   return -2.0*pow(z, 2.0)+pow(x, 2.0)+pow(y, 2.0);
}

//-----------------------------------------------------------------------------

double spherical_harmonic_function_scaled(double x, double y, double z)
{
    //    double l = 4.0;
    double l = 3.0;
    //    double l = 2.0;

    double evalue = -l * (l + 1.0);

    return spherical_harmonic_function(x, y, z) / evalue;
}

//-----------------------------------------------------------------------------

void solve_laplace_equation(SurfaceMesh &mesh, int face_point)
{
    Eigen::SparseMatrix<double> M, S, S_f;
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

    setup_stiffness_matrix(mesh, S, face_point);
    int nb = 0;
    for (auto v : mesh.vertices())
    {
        //count nr outer vertices
        if (mesh.is_boundary(v))
        {
            nb++;
        }
    }
    Eigen::SparseMatrix<double> G, V, Div, P;

    unsigned ins = 0;
    unsigned out = 0;

    Eigen::VectorXi in(mesh.n_vertices() - nb), bound(nb), x(3);

    x << 0, 1, 2;

    for (auto v : mesh.vertices())
    {
        // save indices of inner and outer vertices
        if (!mesh.is_boundary(v))
        {
            in(ins) = v.idx();
            ins++;
        }
        else
        {
            bound(out) = v.idx();
            out++;
            // right-hand side: fix x coordinate of boundary vertices for the righthandsite
            B(v.idx(), 0) = mesh.position(v)[0];
            B(v.idx(), 1) = mesh.position(v)[1];
            B(v.idx(), 2) = mesh.position(v)[2];
        }
    }

    Eigen::SparseMatrix<double> L_in_in, L_in_b;
    Eigen::MatrixXd b_in, b_out;

    // slice S and b and bring boundary values on the righthandsite to solve only for inner vertices
    igl::slice(S, in, in, L_in_in);
    igl::slice(S, in, bound, L_in_b);
    igl::slice(B, in, x, b_in);
    igl::slice(B, bound, x, b_out);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L_in_in);
    Eigen::MatrixXd X = solver.solve(b_in - L_in_b * b_out);
    double error = 0.0;

    int k = 0;
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "harmonic(): Could not solve linear system\n";
    }
    else
    {
        // copy solution
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                error += pow(X(k, 0) - mesh.position(v)[0], 2) +
                         pow(X(k, 1) - mesh.position(v)[1], 2) +
                         pow(X(k, 2) - mesh.position(v)[2], 2);
                k++;
            }
        }
    }
    std::cout << "RMSE inner verticex positions : " << sqrt(error / (double)k)
              << std::endl;
}
