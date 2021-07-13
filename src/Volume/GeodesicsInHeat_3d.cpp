//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "GeodesicsInHeat_3d.h"
#include "LaplaceConstruction_3D.h"
#include "Diamond_3D.h"
#include "diffgeo_3D.h"
#include <Eigen/Sparse>
#include <float.h>

//=============================================================================

GeodesicsInHeat_3d::GeodesicsInHeat_3d(VolumeMesh &mesh,
                                       unsigned int face_point,
                                       unsigned int cell_point)
    : mesh_(mesh), face_point_(face_point), cell_point_(cell_point)
{
}

//-----------------------------------------------------------------------------

GeodesicsInHeat_3d::~GeodesicsInHeat_3d() {}

//-----------------------------------------------------------------------------

double GeodesicsInHeat_3d::maxEdgeLength(const VolumeMesh &mesh)
{
    double maxLen = 0.;
    double currLen;
    for (auto e : mesh.edges())
    {
        currLen = mesh.length(e);
        if (currLen > maxLen)
        {
            maxLen = currLen;
        }
    }
    return maxLen;
}

double GeodesicsInHeat_3d::meanEdgeLength(const VolumeMesh &mesh)
{
    double meanLen = 0.;
    for (auto e : mesh.edges())
    {
        meanLen += mesh.length(e);
    }
    meanLen /= (double)mesh.n_edges();
    return meanLen;
}
//-----------------------------------------------------------------------------

double GeodesicsInHeat_3d::compute_geodesics(const int vertex,
                                             Eigen::VectorXd &dist,
                                             Eigen::VectorXd &orthodist,
                                             bool max_e_len)
{
    Eigen::SparseMatrix<double> L, M, A, M_bar, Grad, Div, P;
    setup_3D_stiffness_matrix(mesh_, L, face_point_, cell_point_);
    setup_3D_mass_matrix(mesh_, M, face_point_, cell_point_);
    Eigen::SparseMatrix<double> V, Pc, Pf, G;
    auto f_w_prop =
        mesh_.request_face_property<Eigen::VectorXd>("face weights");
    auto c_w_prop =
        mesh_.request_cell_property<Eigen::VectorXd>("cell weights");
    auto f_prop =
        mesh_.request_face_property<VolumeMesh::PointT>("face points");
    auto c_prop =
        mesh_.request_cell_property<VolumeMesh::PointT>("cell points");
    compute_3D_virtual_points_and_weights(mesh_, face_point_, cell_point_);
    setup_3D_cell_face_prolongation_matrix(mesh_, P, Pc, Pf);
    setup_3D_diamond_gradient(mesh_, G, V);
    Grad = G * P;
    Div = P.transpose() * -G.transpose() * V;
    double h;
    if(max_e_len){
        h= pow(maxEdgeLength(mesh_), 2);

    }else{
        h = pow(meanEdgeLength(mesh_), 2);
    }

    A = M - h * L;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholL;
    cholL.analyzePattern(L);
    cholL.factorize(L);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholA;
    cholA.analyzePattern(A);
    cholA.factorize(A);

    // diffuse heat
    const int N = mesh_.n_vertices();

    auto distances = mesh_.request_vertex_property<Scalar>("v:dist");

    Eigen::SparseVector<double> b(N);
    b.coeffRef(vertex) = 1.;
    // compute gradients

    Eigen::VectorXd heat = cholA.solve(b);
    Eigen::VectorXd grad = Grad * heat;

    // normalize gradients
    for (int i = 0; i < grad.rows(); i += 3)
    {
        dvec3 &g = *reinterpret_cast<dvec3 *>(&grad[i]);
        double n = norm(g);
        if (n > DBL_MIN)
            g /= n;
    }

    dist = cholL.solve(Div * (-grad));

    orthodist.resize(dist.size());

    double mi = dist.minCoeff();
    for (int i = 0; i < dist.rows(); ++i)
        dist[i] -= mi;
    int k = 0;

    OpenVolumeMesh::VertexHandle v0 = OpenVolumeMesh::VertexHandle(vertex);
    double rms = 0.0;

    for (auto v : mesh_.vertices())
    {
        distances[v] = dist(v.idx());

        //        std::cout << "Result: " << dist(v.idx()) << " Euklid (v-v0): "
        //                  << (mesh_.vertex(v0) - mesh_.vertex(v)).norm() << std::endl;
        orthodist(k) = (mesh_.vertex(v0) - mesh_.vertex(v)).norm();
        rms += (dist(k) - orthodist(k)) * (dist(k) - orthodist(k));

        k++;
    }
    rms /= mesh_.n_vertices();
    rms = sqrt(rms);
    distance_to_texture_coordinates();
    std::cout << "Diamond rmse: " << rms << std::endl;
    std::cout << "-------------------" << std::endl;
    return rms;
}

void GeodesicsInHeat_3d::distance_to_texture_coordinates() const
{

    auto distances = mesh_.request_vertex_property<Scalar>("v:dist");
    assert(distances);
    // scale with boundingbox size for comparable geodesic rings
    Scalar b_sphere_diameter = 2.0*mesh_.get_bounding_sphere().second;
    auto tex = mesh_.request_vertex_property<vec2>("texture_coordinates");
    for (auto v : mesh_.vertices())
    {
        if (distances[v] <= FLT_MAX)
        {
            // std::cout << distances[v] / b_sphere_diameter << std::endl;
            tex[v] = vec2(distances[v] / b_sphere_diameter, 0.0);
        }
        else
        {
            tex[v] = vec2(1.0, 0.0);
        }
    }

    // remove per-halfedge texture coordinates
//    auto htex = mesh_.get_halfedge_property<TexCoord>("h:tex");
//    if (htex)
//        mesh_.remove_halfedge_property(htex);
}