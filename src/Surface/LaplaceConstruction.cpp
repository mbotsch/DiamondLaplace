//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "LaplaceConstruction.h"
#include "diffgeo.h"
#include "Diamond_2D.h"

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

using namespace std;

enum InsertedPoint
{
    Centroid = 0,
    AreaMinimizer = 2,
};

//----------------------------------------------------------------------------------

void setup_stiffness_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S,
                            int minpoint)
{

    if (!mesh.has_face_property("f:point"))
    {
        mesh.add_face_property<Point>("f:point");
        mesh.add_face_property<Eigen::VectorXd>("f:weights");
    }
    setup_face_point_properties(mesh, minpoint);
    Eigen::SparseMatrix<double> G, D, G2, D2, Gra, Div, P;
    setup_prolongation_matrix(mesh, P);
    compute_primal_points(mesh, minpoint);
    setup_diamond_gradient_divergence_intrinsic(mesh, G, D);
    Gra = G * P;
    Div = P.transpose() * D;
    S = Div * Gra;
}

//----------------------------------------------------------------------------------

void setup_mass_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M,
                       int minpoint)
{
    if (!mesh.has_face_property("f:point"))
    {
        mesh.add_face_property<Point>("f:point");
        mesh.add_face_property<Eigen::VectorXd>("f:weights");
    }
    setup_face_point_properties(mesh, minpoint);
    Eigen::SparseMatrix<double> M_, P;
    setup_prolongation_matrix(mesh, P);
    setup_diamond_mass_matrix(mesh, M_);
    M = P.transpose() * M_ * P;
}
