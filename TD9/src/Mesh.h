#ifndef MESH_H
#define MESH_H
#include "Eigen/Core"
#include <vector>

//forward declarations
#include "MeshParts.h"
using namespace MeshParts;


/**
 * @class Mesh
 * @brief Represents a mesh consisting of vertices, faces, and other related properties.
 * 
 * The Mesh class provides functionality to compute various properties of the mesh, such as
 * circumcenters, normals, tangent spaces, and more. It also stores information about the
 * primal and dual vertices, faces, edges, and half-edges of the mesh.
 */
class Mesh
{
private:
    /* data */
    Eigen::MatrixXd compute_circumcenter_triangle(const Eigen::MatrixXd, const Eigen::MatrixXd, const Eigen::MatrixXd);
    void initialize_complex();
    void compute_half_edges();

protected:
    

    //helper function to facilitate the computations
    Eigen::MatrixXd toMatrix(Eigen::Vector3d v){Eigen::MatrixXd M = (Eigen::MatrixXd(1,3)<<v[0],v[1],v[2]).finished();return M;}
    //helper function to facilitate the computations
    Eigen::Vector3d toVector(Eigen::MatrixXd M){Eigen::Vector3d v = Eigen::Vector3d( M(0,0),M(0,1),M(0,2));return v;}

    

public:
    
    Eigen::MatrixXd V; // primal vertices
    Eigen::MatrixXi F; // primal triangular faces
    
    int euler_characteristic;
    int boundaries = -1;


    Eigen::MatrixXd N_faces; // normals of the faces
    Eigen::MatrixXd V_dual;  // dual vertices
    std::vector<std::vector<int>> F_dual; // dual polytopal faces
    
    std::vector<VertexPtr> primal_vertices;
    std::vector<PrimalFacePtr> primal_faces;
    std::vector<HalfEdgePtr> hedges;

    Mesh(/* args */);
    Mesh(Eigen::MatrixXd, Eigen::MatrixXi);

    Eigen::VectorXd compute_gaussian_curvature();
    void compute_voronoi_area();
    void vertexDegreeStatistics();

    Eigen::MatrixXd compute_vertex_normals();
    Eigen::MatrixXd compute_vertex_normals_hed();

    Eigen::MatrixXd compute_face_normals();
    Eigen::MatrixXd compute_face_normals_hed();

    void count_boundaries();
    //------------------------------------------

    ~Mesh();
};

#endif