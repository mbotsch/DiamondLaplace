//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "Eigenmodes.h"
#include <cmath>
#include "LaplaceConstruction_3D.h"
#include <igl/slice.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/Util/GEigsMode.h>
#include <fstream>

using namespace Spectra;

//=============================================================================

enum VolumePoints
{
    Quadratic_Volume_ = 0,
    Cell_Centroid_ = 1
};

double solve_eigenvalue_problem(VolumeMesh &mesh, Eigen::VectorXd &evalues,
                                int face_point, int cell_point,
                                std::string meshname)
{

    std::string filename;

    if (cell_point == Quadratic_Volume_)
    {
        filename = "eigenmodes_Diamond-VolMinimizer-" + meshname + ".csv";
    }
    else
    {
        filename = "eigenmodes_Diamond-Centroid-" + meshname + ".csv";
    }

    std::ofstream ev_file(filename);
    ev_file << "computed,analytic,offset" << std::endl;

    Eigen::SparseMatrix<double> M, S;
    Eigen::VectorXd b(mesh.n_vertices());

    setup_3D_stiffness_matrix(mesh, S, face_point, cell_point);
    setup_3D_mass_matrix(mesh, M, face_point, cell_point);

    // collect indices of inner vertices
    std::vector<int> indices;
    for (auto v : mesh.vertices())
    {
        if (!mesh.is_boundary(v))
        {
            indices.push_back(v.idx());
        }
    }

    Eigen::VectorXi in(indices.size());

    //Rewrite indices to Eigen::Vector
    for (unsigned int i = 0; i < indices.size(); ++i)
    {
        in(i) = indices[i];
    }

    Eigen::SparseMatrix<double> S_in_in, M_in_in;

    //slice matrices so that only rows and cols for inner vertices remain

    igl::slice(S, in, in, S_in_in);
    igl::slice(M, in, in, M_in_in);

    // Construct matrix operation object using the wrapper class SparseGenMatProd
    SparseSymMatProd<double> op(S_in_in);
    SparseCholesky<double> Bop(M_in_in);
    int num_eval = 34;
    int converge_speed = 5 * num_eval;

    // Construct generalized eigen solver object, requesting the smallest generalized eigenvalues
    SymGEigsSolver<
        SparseSymMatProd<double>,
        SparseCholesky<double>, 
        GEigsMode::Cholesky>
        geigs(op, Bop, num_eval, converge_speed);

    // Initialize and compute
    geigs.init();
    geigs.compute(SortRule::SmallestMagn);

    // Retrieve results
    Eigen::VectorXd evectors, analytic;
    if (geigs.info() == CompInfo::Successful)
    {
        evalues = geigs.eigenvalues();
        evectors = geigs.eigenvectors();
    }
    analytic_eigenvalues_unitBall(analytic, num_eval);
    double error = 0.0;
    for (int i = 0; i < evalues.size(); i++)
    {
        ev_file << -evalues(i) << "," << analytic(i) << ","
                << -evalues(i) - analytic(i) << std::endl;
        std::cout << "Computed evalue: " << -evalues(i)
                  << " analytical Bessel: " << analytic(i) << std::endl;

        error += pow(-evalues(i) - analytic(i), 2);
    }
    error = sqrt(error / (double)evalues.size());
    std::cout << "Root mean squared error: " << error << std::endl;
    return error;
}

void analytic_eigenvalues_unitBall(Eigen::VectorXd &eval, int n)
{
    //n so far only <= 34
    if (n > 34)
    {
        std::cout << "n must be lower than 34! \n";
    }
    eval.resize(n);
    for (int i = 0; i < n; i++)
    {
        if (i == 0)
        {
            //Bessel j zero (1/2,1)
            eval(i) = 9.8696;
        }
        else if (i > 0 && i < 4)
        {
            //Bessel j zero (3/2,1)

            eval(i) = 20.1907;
        }
        else if (i >= 4 && i < 9)
        {
            //Bessel j zero (5/2,1)

            eval(i) = 33.2175;
        }
        else if (i == 9)
        {
            //Bessel j zero (1/2,2)

            eval(i) = 39.4784;
        }
        else if (i > 9 && i < 17)
        {
            //Bessel j zero (7/2,1)
            eval(i) = 48.8312;
        }
        else if (i >= 17 && i < 20)
        {
            //Bessel j zero (3/2,2)
            eval(i) = 59.6795;
        }
        else if (i >= 20 && i < 29)
        {
            //Bessel j zero (9/2,1)
            eval(i) = 66.9543;
        }
        else if (i >= 29 && i < 34)
        {
            //Bessel j zero (5/2,2)
            eval(i) = 82.7192;
        }
    }
}
