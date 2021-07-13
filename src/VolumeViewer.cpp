//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include <imgui.h>
#include <VolumeSubdivision.h>
#include "VolumeViewer.h"
#include "Volume/Eigenmodes.h"
#include "Volume/Franke_PoissonSystem_3D.h"
#include "Volume/GeodesicsInHeat_3d.h"

//=============================================================================

enum VolumePoints
{
    Quadratic_Volume_ = 0,
    Cell_Centroid_ = 1
};

enum AreaPoints
{
    Quadratic_Areas_ = 0,
    Face_Centroid = 1
};

bool VolumeViewer::load_mesh(const char *filename)
{
    if (VolumeMeshViewer::load_mesh(filename))
    {
        return true;
    }
    return false;
}

void VolumeViewer::keyboard(int key, int code, int action, int mod)
{
    if (action != GLFW_PRESS) // only react on key press events
        return;

    switch (key)
    {

        case GLFW_KEY_B: {
            auto mesh = static_cast<VolumeMesh>(mesh_);

            mat3 rotation(0.0);
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    rotation(i, j) = modelview_matrix_(i, j);
                }
            }

            {
                auto box = mesh.get_bounding_sphere();

                Vec3f cur = Vec3f(0.0);
                Vec3f max = Vec3f(std::numeric_limits<float>::min());
                Vec3f min = Vec3f(std::numeric_limits<float>::max());
                for (auto v : mesh.vertices())
                {
                    auto vec = mesh.vertex(v);
                    vec -= box.first;
                    auto pmp_vec = pmp::vec3(vec[0], vec[1], vec[2]);
                    pmp_vec = rotation * pmp_vec;
                    vec = Vec3f(pmp_vec[0], pmp_vec[1], pmp_vec[2]);

                    cur = vec;
                    for (int i = 0; i < 3; ++i)
                    {
                        if (cur[i] > max[i])
                        {
                            max[i] = cur[i];
                        }
                        if (cur[i] < min[i])
                        {
                            min[i] = cur[i];
                        }
                    }
                    mesh.set_vertex(v, vec);
                }

                Vec3f diag = max - min;
                float maxD = *std::max_element(diag.data(), diag.data() + 3);

                for (auto v : mesh.vertices())
                {
                    auto vec = mesh.vertex(v);
                    vec /= maxD;
                    mesh.set_vertex(v, vec);
                }
            }

            std::string name = filename_;
            int last_slash = name.find_last_of('/');
            int last_dot = name.find_last_of('.');
            name = name.substr(last_slash + 1, last_dot - last_slash - 1);

            file_manager_.writeFile(name + ".ovm", mesh);
            std::cout << "Created " << name << ".ovm" << std::endl;
            break;
        }

        default: {
            VolumeMeshViewer::keyboard(key, code, action, mod);
            break;
        }
    }
}

void VolumeViewer::process_imgui()
{
    VolumeMeshViewer::process_imgui();

    static int face_point = 0;
    ImGui::RadioButton("Area minimizer", &face_point, Quadratic_Areas_);
    ImGui::RadioButton("Face Centroid", &face_point, Face_Centroid);
    face_point_ = (unsigned)face_point;
    ImGui::Spacing();
    ImGui::Spacing();

    static int cell_point = 0;
    ImGui::RadioButton("Volume minimizer", &cell_point, Quadratic_Volume_);
    ImGui::RadioButton("Cell Centroid", &cell_point, Cell_Centroid_);
    cell_point_ = (unsigned)cell_point;

    ImGui::Spacing();

    static int ts = 1;
    ImGui::RadioButton("Mean edge length", &ts, 0);
    ImGui::RadioButton("Max edge length", &ts, 1);
    ts_ = ts;

    if (ImGui::Button("RMSE Eigenvalues Sphere"))
    {
        Eigen::VectorXd evalues;
        solve_eigenvalue_problem(mesh_, evalues, face_point_, cell_point_);
    }
    if (ImGui::Button("RMSE Franke Poisson System"))
    {
        Eigen::VectorXd results;
        solve_franke_poisson(mesh_, results, face_point_, cell_point_);
    }
    if (ImGui::Button("Linear Precision"))
    {
        Eigen::VectorXd results;
        solve_laplace_equation(mesh_, face_point_, cell_point_);
    }
    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Polyhedra!"))
    {
        if (ImGui::Button("tetrahedra"))
        {
            VolumeSubdivision(mesh_).tetrahedra();
            update_mesh();
        }

        if (ImGui::Button("irregular pyrmaids"))
        {
            VolumeSubdivision(mesh_).irregular_mesh(5);
            update_mesh();
        }

        if (ImGui::Button("full truncation"))
        {
            VolumeSubdivision(mesh_).full_truncation();
            update_mesh();
        }

        if (ImGui::Button("Quads"))
        {
            VolumeSubdivision(mesh_).quads();
            update_mesh();
        }

        if (ImGui::Button("Linear Subdivision"))
        {
            VolumeSubdivision(mesh_).linear_subdivision();
            update_mesh();
        }
    }
    if (ImGui::CollapsingHeader("Geodesics in Heat",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Compute Distances Vertex 0"))
        {
            GeodesicsInHeat_3d heat(mesh_, face_point_, cell_point_);
            Eigen::VectorXd dist, geodist;
            heat.compute_geodesics(0, dist, geodist, ts_);
            update_mesh();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");
        }
    }
}

void VolumeViewer::mouse(int button, int action, int mods)
{

    if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_MIDDLE &&
        mods == GLFW_MOD_SHIFT)
    {
        double x, y;
        cursor_pos(x, y);
        VHandle v = pick_vertex(x, y);
        std::cout << "Vertex Idx: " << v.idx() << std::endl;
        if (v.idx() != -1)
        {
            GeodesicsInHeat_3d heat(mesh_, face_point_, cell_point_);
            Eigen::VectorXd dist, geodist;
            heat.compute_geodesics(v.idx(), dist, geodist, ts_);
            //            heat.distance_to_texture_coordinates();
            update_mesh();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");
        }
    }
    else
    {
        VolumeMeshViewer::mouse(button, action, mods);
    }
}
