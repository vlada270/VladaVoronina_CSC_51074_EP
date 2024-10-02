#include "Mesh.h"
#include <iostream>
#include "HalfedgeDS.h"
#include "HalfedgeBuilder.h"
#include "igl/edges.h"
#include <algorithm>
#include <chrono>
#include <memory>
#include "igl/gaussian_curvature.h"
#include "igl/massmatrix.h"


using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Vector3d;
using Eigen::VectorXd;

Mesh::Mesh(/* args */)
{
}

/**
 * @brief Constructs a Mesh object with the given vertex and face matrices.
 * 
 * @param V The vertex matrix representing the coordinates of the vertices.
 * @param F The face matrix representing the vertex indices of each face.
 */
Mesh::Mesh(MatrixXd V, MatrixXi F)
{
    this->V = V;
    this->F = F;

    //TODO: Implement the constructor

    


    
    std::cout << "Calling constructor of the mesh" << std::endl;
    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();
    

    initialize_complex();
    auto stop0 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop0 - start;
    std::cout << "Elapsed time for complex initialization: " << elapsed.count() << " seconds" << std::endl;

    compute_half_edges();
    auto stop1 = std::chrono::high_resolution_clock::now();
    elapsed = stop1 - stop0;
    std::cout << "Elapsed time for hedge initialization: " << elapsed.count() << " seconds" << std::endl;

    
}

Mesh::~Mesh()
{
}

/**
 * Computes the circumcenter of a triangle given its three vertices.
 *
 * @param v0 The coordinates of the first vertex.
 * @param v1 The coordinates of the second vertex.
 * @param v2 The coordinates of the third vertex.
 * @return The coordinates of the circumcenter as a 1x3 matrix.
 */
MatrixXd Mesh::compute_circumcenter_triangle(const MatrixXd v0, const MatrixXd v1, const MatrixXd v2)
{
    // TODO
}

/**
 * @brief Initializes the complex by creating primal faces and vertices.
 * 
 * This function calculates the barycenter, circumcenter, normal, and tangent space basis for each primal face.
 * It also initializes the primal vertices with their positions.
 */
void Mesh::initialize_complex()
{
    //TODO: Implement the initialization of the complex

}

/**
 * Computes the half edges of the mesh.
 * 
 * This function loops over the ring for each vertex and creates the half edges
 * based on the vertex and face information. It also sets up the necessary
 * connections between the half edges, vertices, edges, and faces.
 * 
 * @return void
 */
void Mesh::compute_half_edges()
{

    //TODO


}

//------------------------------------


void Mesh::vertexDegreeStatistics() {
		
    std::cout << "Computing vertex degree distribution " << std::endl;
    auto start = std::chrono::high_resolution_clock::now(); // for measuring time performances
    int *vDS = new int[V.rows()];
    //TODO: Implement the vertex degree statistics
    
    for (int i = 3; i < V.rows(); i++)
    	if (vDS[i] != 0)
    		std::cout << "number of degrees = " << i << " ; number of occurences = " << vDS[i] << std::endl;
    std::cout << "Zero for the other vertex degrees between d=3 and n-1" << std::endl;
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
}

Eigen::VectorXd Mesh::compute_gaussian_curvature(){
    //TODO
}

void Mesh::compute_voronoi_area(){
    //TODO
}

void Mesh::count_boundaries(){
    //TODO
}

MatrixXd Mesh::compute_vertex_normals(){
    //TODO
}

MatrixXd Mesh::compute_vertex_normals_hed(){
    //TODO
}

MatrixXd Mesh::compute_face_normals(){
    //TODO
}

MatrixXd Mesh::compute_face_normals_hed(){
    //TODO
}