//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include "../SurfaceViewer.h"
#include <pmp/SurfaceMesh.h>
#include <Eigen/Sparse>
#include <fstream>
#include <pmp/algorithms/SurfaceNormals.h>
//=============================================================================

class GeodesicsInHeat
{
public:
    GeodesicsInHeat(pmp::SurfaceMesh& mesh, int min_point, bool geodist,
                    bool euklid, bool mean_edge_ = false);
    ~GeodesicsInHeat();

    void getDistance(const int vertex, Eigen::VectorXd& dist,
                     Eigen::VectorXd& orthodist);

    void distance_to_texture_coordinates() const;

    void compute_geodesics();

private:
    SurfaceMesh& mesh_;

    unsigned int min_point_;

    Eigen::MatrixXd pos;

    Eigen::SparseMatrix<double> divOperator, gradOperator;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholL, cholA;

    double averageEdgeLength(const pmp::SurfaceMesh& mesh);

    double maxEdgeLength(const pmp::SurfaceMesh& mesh);

    double great_circle_distance(Vertex v, Vertex vv, double r = 1.0);

    double haversine_distance(Vertex v, Vertex vv, double r = 1.0);

    double vincenty_distance(Vertex v, Vertex vv, double r = 1.0);

    bool geodist_sphere_, geodist_cube_;

    bool mean_edge_;
};
