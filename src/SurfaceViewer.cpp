//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "SurfaceViewer.h"
#include "Surface/Parameterization.h"
#include "Surface/Smoothing.h"
#include "Surface/diffgeo.h"
#include "Surface/Curvature.h"
#include "Surface/SpectralProcessing.h"
#include "Surface/GeodesicsInHeat.h"
#include "Surface/Diamond_2D.h"
#include "Surface/Poisson_System.h"

#include <pmp/algorithms/SurfaceTriangulation.h>
#include <pmp/algorithms/SurfaceSubdivision.h>
#include <pmp/algorithms/SurfaceNormals.h>
#include <pmp/algorithms/SurfaceGeodesic.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/Timer.h>

#include <imgui.h>
#include <unsupported/Eigen/SparseExtra>
#include <random>

using namespace pmp;

//=============================================================================

enum InsertedPoint
{
    Centroid_ = 0,
    AbsAreaMinimizer = 1,
    AreaMinimizer = 2,
    Triangle_Circumcenter = 3
};

void Viewer::keyboard(int key, int scancode, int action, int mods)
{
    if (action != GLFW_PRESS) // only react on key press events
        return;

    switch (key)
    {

        case GLFW_KEY_L: {
            std::string line;
            std::ifstream myfile;
            myfile.open("modelview_matrix.txt");
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    getline(myfile, line);
                    modelview_matrix_(i, j) = std::stof(line);
                }
            }
            myfile.close();
            break;
        }

        case GLFW_KEY_S: {
            std::ofstream myfile;
            myfile.open("modelview_matrix.txt");
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    myfile << modelview_matrix_(i, j) << "\n";
                    std::cout << modelview_matrix_(i, j) << "\n";
                }
            }
            myfile.close();
            break;
        }

        case GLFW_KEY_B: {
            auto mesh2 = static_cast<SurfaceMesh>(mesh_);
            auto points = mesh2.get_vertex_property<Point>("v:point");

            mat3 rotation(0.0);
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    rotation(i, j) = modelview_matrix_(i, j);
                }
            }

            {
                BoundingBox box;

                for (auto v : mesh2.vertices())
                    box += points[v];

                auto diag = (box.max() - box.min());
                float maxD = *std::max_element(diag.data(), diag.data() + 3);

                for (auto v : mesh2.vertices())
                {
                    points[v] -= box.center();
                    points[v] = rotation * points[v];
                    points[v] /= maxD;
                }
            }

            std::string name = filename_;
            int last_slash = name.find_last_of('/');
            int last_dot = name.find_last_of('.');
            name = name.substr(last_slash + 1, last_dot - last_slash - 1);

            auto vtex = mesh2.get_vertex_property<TexCoord>("v:tex");
            auto htex = mesh2.get_halfedge_property<TexCoord>("h:tex");

            if (vtex)
            {
                if (!htex)
                {
                    htex = mesh2.add_halfedge_property<TexCoord>("h:tex");
                }
                for (SurfaceMesh::HalfedgeIterator hit =
                         mesh2.halfedges_begin();
                     hit != mesh2.halfedges_end(); ++hit)
                {
                    htex[*hit] = vtex[mesh2.to_vertex(*hit)];
                }
            }

            SurfaceNormals::compute_vertex_normals(mesh2);
            mesh2.write(name + ".obj");
            std::cout << "Created " << name << ".obj" << std::endl;
            break;
        }
        default: {
            MeshViewer::keyboard(key, scancode, action, mods);
            break;
        }
    }
}

//----------------------------------------------------------------------------

void Viewer::process_imgui()
{
    // add standard mesh info stuff
    pmp::MeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Settings", ImGuiTreeNodeFlags_DefaultOpen))
    {

        ImGui::Text("Choose your minimizing Point ");

        ImGui::Spacing();

        static int min_point = 2;

        ImGui::RadioButton("Centroid", &min_point, 0);
        ImGui::RadioButton("Area Minimizer", &min_point, 2);

        ImGui::Text("Choose your diffusion timestep ");

        ImGui::Spacing();

        static int ts = 0;

        ImGui::RadioButton("Mean edge length", &ts, 0);
        ImGui::RadioButton("Max edge length", &ts, 1);

        min_point_ = min_point;
        if (ts == 0)
        {
            time_step_ = true;
        }
        else
        {
            time_step_ = false;
        }
    }

    ImGui::Spacing();

    if (ImGui::Button("Linear Precision"))
    {
        solve_laplace_equation(mesh_, min_point_);
    }

    ImGui::Spacing();
    ImGui::Spacing();
    // turn mesh into non-triangles
    if (ImGui::CollapsingHeader("Polygons!"))
    {
        // Catmull-Clark subdivision
        if (ImGui::Button("Catmull-Clark"))
        {
            pmp::SurfaceSubdivision subdiv(mesh_);
            subdiv.catmull_clark();
            update_mesh();
        }
        // loop subdivision
        if (ImGui::Button("Loop"))
        {
            SurfaceSubdivision(mesh_).loop();
            update_mesh();
        }
        if (ImGui::Button("insert virtual points"))
        {
            insert_points(min_point_);
        }
        if (ImGui::Button("Centroids"))
        {
            Centroid();
        }

        // dualize the mesh
        if (ImGui::Button("Dualize mesh"))
        {
            dualize();
        }

        // close holes by polygons
        if (ImGui::Button("Close holes"))
        {
            close_holes();
        }

        if (ImGui::Button("Triangulate mesh (min area)"))
        {
            SurfaceTriangulation tesselator(mesh_);
            tesselator.triangulate(SurfaceTriangulation::Objective::MIN_AREA);
            update_mesh();
        }

        if (ImGui::Button("Triangulate mesh (max angle)"))
        {
            SurfaceTriangulation tesselator(mesh_);
            tesselator.triangulate(SurfaceTriangulation::Objective::MAX_ANGLE);
            update_mesh();
        }
        if (ImGui::Button("Kugelize"))
        {
            for (auto v : mesh_.vertices())
                mesh_.position(v) = normalize(mesh_.position(v));
            update_mesh();
        }
        if (ImGui::Button("Add noise"))
        {
            for (auto v : mesh_.vertices())
            {
                Point n = SurfaceNormals::compute_vertex_normal(mesh_, v);
                Scalar r = 2.0 * static_cast<float>(rand()) /
                               static_cast<float>(RAND_MAX) -
                           1.0;
                Scalar h = r * 0.01 * radius_;
                mesh_.position(v) += h * n;
            }
            update_mesh();
        }
        static bool normal_, fixed_;

        ImGui::Checkbox("noise in normal direction?", &normal_);
        ImGui::Checkbox("Constant normal noise?", &fixed_);
        if (ImGui::Button("Noise"))
        {
            // create random generator
            // upper and lower bounds are proportional to bounding box and inverse of smoothness
            auto l = mesh_.bounds().size();
            double upper_bound = l / 1000.0;
            double lower_bound = -upper_bound;
            std::uniform_real_distribution<double> unif(lower_bound,
                                                        upper_bound);
            std::default_random_engine re;

            re.seed(42); // fixed seed

            double rand = 0;
            if (fixed_)
            {
                rand = unif(re);
            }

            for (auto v : mesh_.vertices())
            {
                auto n = pmp::SurfaceNormals::compute_vertex_normal(mesh_, v);
                if (normal_)
                {
                    if (fixed_)
                    {
                        mesh_.position(v) -= n * rand;
                    }
                    else
                    {
                        mesh_.position(v) -= n * unif(re);
                    }
                }
                else
                {
                    mesh_.position(v) +=
                        pmp::Point(unif(re), unif(re), unif(re));
                }
            }

            update_mesh();
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();

    // discrete harmonic parameterization
    if (ImGui::CollapsingHeader("Parametrization"))
    {
        if (ImGui::Button("Discrete Harmonic"))
        {
            Parameterization(mesh_).harmonic(min_point_);
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");
            update_mesh();
        }
    }
    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Poisson System",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::PushItemWidth(100);
        static int function = 1;
        ImGui::RadioButton("Franke 2D (planar)", &function, 0);
        ImGui::RadioButton("Spherical Harmonics", &function, 1);

        ImGui::PopItemWidth();

        if (ImGui::Button("Solve!"))
        {
            solve_poisson_system(mesh_, min_point_, function);
        }
    }
    ImGui::Spacing();
    ImGui::Spacing();

    // implicit smoothing
    if (ImGui::CollapsingHeader("Smoothing", ImGuiTreeNodeFlags_DefaultOpen))
    {
        static float timestep = 0.1;
        float lb = 0.001;
        float ub = 1.0;
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("TimeStep", &timestep, lb, ub);
        ImGui::PopItemWidth();

        if (ImGui::Button("Implicit Smoothing"))
        {
            close_holes();

            Scalar dt = timestep;
            smooth_.implicit_smoothing(dt, min_point_);
            update_mesh();
            BoundingBox bb = mesh_.bounds();
            set_scene((vec3)bb.center(), 0.5 * bb.size());
            open_holes();
        }
        if (ImGui::Button("20 * Implicit Smoothing"))
        {
            close_holes();

            Scalar dt = timestep;
            for (int i = 0; i < 20; i++)
            {
                smooth_.implicit_smoothing(dt, min_point_);
                update_mesh();
            }
            BoundingBox bb = mesh_.bounds();
            set_scene((vec3)bb.center(), 0.5 * bb.size());
            open_holes();
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::Button("Spherical Harmonics"))
    {
        SpectralProcessing analyzer_(mesh_);
        analyzer_.update_eigenvectors(min_point_);
        analyzer_.analyze_eigenvectors();
    }

    // curvature visualization
    if (ImGui::CollapsingHeader("Curvature", ImGuiTreeNodeFlags_DefaultOpen))
    {
        static bool curvature_sphere_ = false;
        ImGui::Checkbox("Compare to unit sphere curvatures",
                        &curvature_sphere_);

        if (ImGui::Button("Mean Curvature"))
        {
            Curvature analyzer(mesh_, curvature_sphere_);
            analyzer.visualize_curvature(min_point_);
            mesh_.use_cold_warm_texture();
            update_mesh();
            set_draw_mode("Texture");
        }
    }
    if (ImGui::CollapsingHeader("Geodesics in Heat",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        static bool geodesic_sphere_ = false;
        static bool geodesic_cube_ = false;
        ImGui::Checkbox("Compare distances to arc lengths", &geodesic_sphere_);
        ImGui::Checkbox("Compare to euclidean distances", &geodesic_cube_);
        compare_sphere = geodesic_sphere_;
        compare_cube = geodesic_cube_;
        if (ImGui::Button("Compute Distances Vertex 0"))
        {
            GeodesicsInHeat heat(mesh_, min_point_, compare_sphere,
                                 compare_cube, time_step_);
            Eigen::VectorXd dist, geodist;
            heat.compute_geodesics();
            heat.getDistance(0, dist, geodist);

            update_mesh();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");

            mesh_.use_checkerboard_texture();
            update_mesh();
            set_draw_mode("Texture");
        }
    }
}

//----------------------------------------------------------------------------

void Viewer::draw(const std::string &draw_mode)
{
    // normal mesh draw
    mesh_.draw(projection_matrix_, modelview_matrix_, draw_mode);
}

//----------------------------------------------------------------------------

void Viewer::dualize()
{
    SurfaceMeshGL dual;

    auto fvertex = mesh_.add_face_property<Vertex>("f:vertex");
    for (auto f : mesh_.faces())
    {
        fvertex[f] = dual.add_vertex(centroid(mesh_, f));
    }

    for (auto v : mesh_.vertices())
    {
        if (!mesh_.is_boundary(v))
        {
            std::vector<Vertex> vertices;
            for (auto f : mesh_.faces(v))
                vertices.push_back(fvertex[f]);
            dual.add_face(vertices);
        }
    }

    mesh_ = dual;
    update_mesh();
}

//----------------------------------------------------------------------------

void Viewer::update_mesh()
{
    // re-compute face and vertex normals
    mesh_.update_opengl_buffers();
}

//----------------------------------------------------------------------------

void Viewer::close_holes()
{
    bool finished = false;
    std::vector<Face> holes;
    while (!finished)
    {
        finished = true;

        // loop through all vertices
        for (auto v : mesh_.vertices())
        {
            // if we find a boundary vertex...
            if (mesh_.is_boundary(v))
            {
                // trace boundary loop
                std::vector<Vertex> vertices;
                vertices.push_back(v);
                for (Halfedge h = mesh_.halfedge(v); mesh_.to_vertex(h) != v;
                     h = mesh_.next_halfedge(h))
                {
                    vertices.push_back(mesh_.to_vertex(h));
                }

                // add boudary loop as polygonal face
                Face f = mesh_.add_face(vertices);
                holes.push_back(f);
                // start over
                finished = false;
                break;
            }
        }
    }
    holes_ = holes;
    update_mesh();
}

//----------------------------------------------------------------------------
void Viewer::insert_points(unsigned int minpoint)
{

    Eigen::MatrixXd poly;
    Eigen::VectorXd w;

    for (Face f : mesh_.faces())
    {
        const int n = mesh_.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (Vertex v : mesh_.vertices(f))
        {
            for (int h = 0; h < 3; h++)
            {
                poly.row(i)(h) = mesh_.position(v)[h];
            }
            i++;
        }

        // compute weights for the polygon
        Eigen::Vector3d p;
        if (minpoint == Centroid_)
        {
            int val = poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double)val;
        }
        else
        {
            find_area_minimizer_weights(poly, w);
        }
        Eigen::Vector3d point = poly.transpose() * w;
        Vertex ver = mesh_.add_vertex(Point(point(0), point(1), point(2)));
        mesh_.split(f, ver);
    }
    mesh_.garbage_collection();
    update_mesh();
}

//----------------------------------------------------------------------------

void Viewer::Centroid()
{
    for (Face f : mesh_.faces())
    {
        Vertex hidx = mesh_.add_vertex(centroid(mesh_, f));
        mesh_.split(f, hidx);
    }
    mesh_.garbage_collection();
    update_mesh();
}

void Viewer::open_holes()
{
    for (Face f : holes_)
    {
        mesh_.delete_face(f);
    }
    mesh_.garbage_collection();
    update_mesh();
}

void Viewer::mouse(int button, int action, int mods)
{

    if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_MIDDLE &&
        mods == GLFW_MOD_SHIFT)
    {
        double x, y;
        cursor_pos(x, y);
        Vertex v = pick_vertex(x, y);
        std::cout << "Vertex Idx: " << v.idx() << std::endl;
        if (mesh_.is_valid(v))
        {
            GeodesicsInHeat heat(mesh_, min_point_, compare_sphere,
                                 compare_cube, time_step_);
            Eigen::VectorXd dist, geodist;

            heat.compute_geodesics();
            heat.getDistance(v.idx(), dist, geodist);
            update_mesh();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");
        }
    }
    else
    {
        MeshViewer::mouse(button, action, mods);
    }
}

void Viewer::load_mesh(const char *filename)
{
    try
    {
        MeshViewer::load_mesh(filename);
        set_draw_mode("Hidden Line");
    }
    catch (const IOException &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

//=============================================================================
