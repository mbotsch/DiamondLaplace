//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "VolumeMeshIO.h"
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <clocale>

//=============================================================================

typedef OpenVolumeMesh::FaceHandle FHandle;
typedef std::vector<VHandle> FaceTuple;
typedef std::map<FaceTuple, FHandle> FaceMap;

//=============================================================================

bool VolumeMeshIO::read(VolumeMesh &mesh)
{
    std::setlocale(LC_NUMERIC, "C");

    // clear mesh before reading from file what does it to a volume mesh
    //mesh.clear();
    for (auto v : mesh.vertices())
    {
        mesh.delete_vertex(v);
    }
    mesh.collect_garbage();

    // extract file extension
    std::string::size_type dot(filename_.rfind('.'));
    if (dot == std::string::npos)
        return false;
    std::string ext = filename_.substr(dot + 1, filename_.length() - dot - 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    // extension determines reader
    if (ext == "ovm")
    {
        return read_ovm(mesh);
    }
    else if (ext == "mesh")
    {
        std::cout << "test" << std::endl;
        return read_mesh(mesh);
    }

    return false;
}

//-----------------------------------------------------------------------------

bool VolumeMeshIO::write(const VolumeMesh &mesh)
{
    // extract file extension
    std::string::size_type dot(filename_.rfind('.'));
    if (dot == std::string::npos)
        return false;
    std::string ext = filename_.substr(dot + 1, filename_.length() - dot - 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);

    // extension determines reader
    if (ext == "ovm")
    {
        return write_ovm(mesh);
    }
    else if (ext == "mesh")
    {
        return write_mesh(mesh);
    }

    return false;
}

//-----------------------------------------------------------------------------

bool VolumeMeshIO::read_ovm(VolumeMesh &mesh)
{
    OpenVolumeMesh::IO::FileManager file_manager;
    return file_manager.readFile(filename_, mesh);
}

//-----------------------------------------------------------------------------

bool VolumeMeshIO::write_ovm(const VolumeMesh &mesh)
{
    OpenVolumeMesh::IO::FileManager file_manager;
    return file_manager.writeFile(filename_, mesh);
}

//-----------------------------------------------------------------------------

void read_mesh_vertices(std::istream &_istream, VolumeMesh &_mesh)
{
    int n_vertices;
    _istream >> n_vertices;

    std::string coord;
    Vec3f vec;
    float x;

    for (int i = 0; i < n_vertices; i++)
    {
        _istream >> vec[0] >> vec[1] >> vec[2] >> x;

        _mesh.add_vertex(vec);
    }
}

//-----------------------------------------------------------------------------

void read_mesh_tetrahedra(std::istream &_istream, VolumeMesh &_mesh)
{
    int n_tetrahedra;
    _istream >> n_tetrahedra;

    std::vector<int> idx;
    idx.resize(4);
    float x;

    std::vector<VHandle> c_vertices_;
    std::vector<FaceTuple> tuples;
    std::vector<HFHandle> cell_halffaces;

    FaceMap face_map;

    for (int i = 0; i < n_tetrahedra; i++)
    {
        idx.clear();
        _istream >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> x;

        c_vertices_.clear();
        c_vertices_.emplace_back(VHandle(idx[0] - 1));
        c_vertices_.emplace_back(VHandle(idx[1] - 1));
        c_vertices_.emplace_back(VHandle(idx[2] - 1));
        c_vertices_.emplace_back(VHandle(idx[3] - 1));

        Vec3f midP(0.0, 0.0, 0.0);
        double valence = 0.0;
        for (auto c_vertex : c_vertices_)
        {
            midP += _mesh.vertex(c_vertex);
            valence += 1.0;
        }
        midP /= valence;

        std::sort(c_vertices_.begin(), c_vertices_.end());

        tuples.clear();
        tuples.emplace_back(
            FaceTuple{c_vertices_[0], c_vertices_[1], c_vertices_[2]});
        tuples.emplace_back(
            FaceTuple{c_vertices_[1], c_vertices_[2], c_vertices_[3]});
        tuples.emplace_back(
            FaceTuple{c_vertices_[0], c_vertices_[2], c_vertices_[3]});
        tuples.emplace_back(
            FaceTuple{c_vertices_[0], c_vertices_[1], c_vertices_[3]});

        cell_halffaces.clear();

        for (const auto &it : tuples)
        {

            // Check if face exists for current tuple
            auto f = face_map.find(it);
            if (f == face_map.end())
            {
                // Face does not exist, create it

                // Find right orientation, s.t. normal
                // points inside the cell

                Vec3f e1 = _mesh.vertex(it[1]) - _mesh.vertex(it[0]);
                Vec3f e2 = _mesh.vertex(it[2]) - _mesh.vertex(it[1]);

                // Get face normal (cross product)
                Vec3f n = (e1 % e2).normalize();

                std::vector<VHandle> v_vec;
                v_vec.push_back(it[0]);
                v_vec.push_back(it[1]);
                v_vec.push_back(it[2]);
                FHandle fh = _mesh.add_face(v_vec);

                // Add face to face map
                face_map[it] = fh;

                // Check whether normal points inside cell
                if (((midP - _mesh.vertex(it[0])) | n) > 0.0)
                {

                    // Normal points inside cell, just add half-face 0
                    // Add corresponding half-face to cell definition
                    cell_halffaces.push_back(_mesh.halfface_handle(fh, 0));
                }
                else
                {

                    // Normal points outside cell, just add half-face 1
                    // Add corresponding half-face to cell definition
                    cell_halffaces.push_back(_mesh.halfface_handle(fh, 1));
                }
            }
            else
            {

                // Face exists, find right orientation
                FHandle fh = f->second;

                std::vector<OpenVolumeMesh::HalfEdgeHandle> hes =
                    _mesh.face(fh).halfedges();

                assert(hes.size() == 3);

                Vec3f e1 = _mesh.vertex(_mesh.halfedge(hes[0]).to_vertex()) -
                           _mesh.vertex(_mesh.halfedge(hes[0]).from_vertex());
                Vec3f e2 = _mesh.vertex(_mesh.halfedge(hes[1]).to_vertex()) -
                           _mesh.vertex(_mesh.halfedge(hes[1]).from_vertex());

                Vec3f n = (e1 % e2).normalize();

                if (((midP -
                      _mesh.vertex(_mesh.halfedge(hes[0]).from_vertex())) |
                     n) > 0.0)
                {
                    // Normal points inside cell
                    cell_halffaces.push_back(_mesh.halfface_handle(fh, 0));
                }
                else
                {
                    // Normal points outisde cell
                    cell_halffaces.push_back(_mesh.halfface_handle(fh, 1));
                }
            }
        }
        _mesh.add_cell(cell_halffaces);
    }
}

void read_mesh_hexahedra(std::istream &_istream, VolumeMesh &_mesh)
{
    int n_hexahedra;
    _istream >> n_hexahedra;

    std::vector<int> idx;
    idx.resize(8);
    float x;

    std::vector<VHandle> c_vertices_;
    std::vector<FaceTuple> tuples;
    std::vector<HFHandle> cell_halffaces;

    FaceMap face_map;

    for (int i = 0; i < n_hexahedra; i++)
    {
        idx.clear();
        _istream >> idx[0] >> idx[1] >> idx[2] >> idx[3] >> idx[4] >> idx[5] >>
            idx[6] >> idx[7] >> x;

        c_vertices_.clear();
        for (int j = 0; j < 8; j++)
        {
            c_vertices_.emplace_back(VHandle(idx[j] - 1));
        }

        Vec3f midP(0.0, 0.0, 0.0);
        double valence = 0.0;
        for (auto c_vertex : c_vertices_)
        {
            midP += _mesh.vertex(c_vertex);
            valence += 1.0;
        }
        midP /= valence;

        int hexahedron_indices_[6][4] = {{0, 1, 2, 3}, {5, 4, 7, 6},
                                         {4, 0, 3, 7}, {1, 5, 6, 2},
                                         {4, 5, 1, 0}, {6, 7, 3, 2}};

        tuples.clear();
        FaceTuple tuple;
        for (auto face : hexahedron_indices_)
        {
            tuple.clear();
            for (size_t k = 0; k < 4; ++k)
            {
                tuple.emplace_back(c_vertices_[face[k]]);
            }
            tuples.emplace_back(tuple);
        }

        cell_halffaces.clear();

        for (const auto &it : tuples)
        {
            std::vector<VHandle> key(it);
            std::sort(key.begin(), key.end());
            auto f = face_map.find(key);
            if (f == face_map.end())
            {

                Vec3f e1 = _mesh.vertex(it[1]) - _mesh.vertex(it[0]);
                Vec3f e2 = _mesh.vertex(it[2]) - _mesh.vertex(it[1]);

                Vec3f n = (e1 % e2).normalize();

                FHandle fh = _mesh.add_face(it);

                face_map[key] = fh;

                if (((midP - _mesh.vertex(it[0])) | n) > 0.0)
                {
                    cell_halffaces.push_back(_mesh.halfface_handle(fh, 0));
                }
                else
                {
                    cell_halffaces.push_back(_mesh.halfface_handle(fh, 1));
                }
            }
            else
            {
                FHandle fh = f->second;

                std::vector<OpenVolumeMesh::HalfEdgeHandle> hes =
                    _mesh.face(fh).halfedges();

                assert(hes.size() == 4);

                Vec3f e1 = _mesh.vertex(_mesh.halfedge(hes[0]).to_vertex()) -
                           _mesh.vertex(_mesh.halfedge(hes[0]).from_vertex());
                Vec3f e2 = _mesh.vertex(_mesh.halfedge(hes[1]).to_vertex()) -
                           _mesh.vertex(_mesh.halfedge(hes[1]).from_vertex());

                Vec3f n = (e1 % e2).normalize();

                if (((midP -
                      _mesh.vertex(_mesh.halfedge(hes[0]).from_vertex())) |
                     n) > 0.0)
                {
                    // Normal points inside cell
                    cell_halffaces.push_back(_mesh.halfface_handle(fh, 0));
                }
                else
                {
                    // Normal points outisde cell
                    cell_halffaces.push_back(_mesh.halfface_handle(fh, 1));
                }
            }
        }
        _mesh.add_cell(cell_halffaces);
    }
}

//-----------------------------------------------------------------------------

bool VolumeMeshIO::read_mesh(VolumeMesh &mesh)
{
    auto iff = std::ifstream(filename_, std::ios::in);

    if (!iff.good())
    {
        std::cerr << "Could not open file " << filename_ << " for reading!"
                  << std::endl;
        return false;
    }

    std::string next;
    float value;
    while (iff.good())
    {
        iff >> next;
        if (next.compare("MeshVersionFormatted") == 0)
        {
            iff >> value;
        }
        else if (next.compare("Dimension") == 0)
        {
            iff >> value;
        }
        else if (next[0] == '#')
        {
            getline(iff, next);
        }
        else if (next.compare("Vertices") == 0)
        {
            read_mesh_vertices(iff, mesh);
        }
        else if (next.compare("Tetrahedra") == 0)
        {
            read_mesh_tetrahedra(iff, mesh);
        }
        else if (next.compare("Hexahedra") == 0)
        {
            read_mesh_hexahedra(iff, mesh);
        }
        else if (next.compare("End") == 0)
        {
            break;
        }
        else
        {
            continue;
        }
    }
    iff.close();

    std::cout << "Converted " << mesh.n_vertices() << " vertices," << std::endl;
    std::cout << "\t  " << mesh.n_edges() << " edges," << std::endl;
    std::cout << "\t  " << mesh.n_faces() << " faces," << std::endl;
    std::cout << "\t  " << mesh.n_cells() << " cells!" << std::endl;

    return true;
}

//-----------------------------------------------------------------------------

bool VolumeMeshIO::write_mesh(const VolumeMesh &mesh)
{
    std::cout << "Currently not supported." << std::endl;
    return false;
}
