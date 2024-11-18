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
#include <unordered_set>
#include <set>


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

    for (int i = 0; i < V.rows(); i++)
    {
        Vertex v = Vertex();
        v.index = i;
        v.position = this->toVector(V.row(i));
        std::shared_ptr<Vertex> vertex = std::make_shared<Vertex>(v);
        primal_vertices.push_back(vertex);
    }

    for (int j = 0; j < F.rows(); j++)
    {
        PrimalFace f = PrimalFace();
        f.index = j;
        for (int k = 0; k < F.cols(); k++) {
            int vertex_index = F(j, k);
            f.indices_vertices.push_back(vertex_index);
        }
        std::shared_ptr<PrimalFace> face = std::make_shared<PrimalFace>(f);
        primal_faces.push_back(face);
    }

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
    // Iterate through each face in the mesh
    for (const auto& face : primal_faces)
    {
        Eigen::Vector3d barycenter(0.0, 0.0, 0.0);
        Eigen::Vector3d normal(0.0, 0.0, 0.0);
        for (int vertex_index : face->indices_vertices) {
            auto vertex = primal_vertices[vertex_index];
            face->addVertex(vertex);
            barycenter += vertex->position;
            vertex->addOneRingFace(face);
        }
        barycenter /= face->indices_vertices.size();
        face->barycenter = barycenter;
        Eigen::Vector3d v0 = primal_vertices[face->indices_vertices[0]]->position;
        Eigen::Vector3d v1 = primal_vertices[face->indices_vertices[1]]->position;
        Eigen::Vector3d v2 = primal_vertices[face->indices_vertices[2]]->position;
        normal = (v1 - v0).cross(v2 - v0).normalized();
        face->normal = normal;

//        face->circumcenter = compute_circumcenter_triangle(this->toMatrix(v0), this->toMatrix(v1), this->toMatrix(v2));
    }

}


// * Computes the half edges of the mesh.
// *
// * This function loops over the ring for each vertex and creates the half edges
// * based on the vertex and face information. It also sets up the necessary
// * connections between the half edges, vertices, edges, and faces.
// *
// * @return void
// */

void Mesh::compute_half_edges() {
    // 1. Setup Half-Edge Builder
    std::unique_ptr<HalfedgeBuilder> builder(new HalfedgeBuilder());
    HalfedgeDS he = builder->createMesh(V.rows(), F);

    //2. Half-Edge Creation and Initialization:
    // After construting HalfedgeDS he, you should create according to the indices of he instances of MeshParts::HalfEdge and store shared pointers to these instances in the std::vector<HalfEdgePtr> this->hedges.
    int num_hedges = F.rows() * 3;
    hedges.resize(num_hedges);

    for (int i = 0; i < num_hedges; i++) {
        hedges[i] = std::make_shared<MeshParts::HalfEdge>(i);
    }

    //3. Half-Edge Connection:
    //Now, that the MeshParts::HalfEdge instances are created, you can use the methods of HalfedgeDS ( read the header file to see what methods are available) to store on each half edge the shared pointers to the start vertex, end vertex, next half edge, flip half edge, boundary etc
    for (int face_idx = 0; face_idx < F.rows(); face_idx++) {
        auto face = primal_faces[face_idx];
        for (int j = 0; j < 3; j++) {
            int hedge_idx = face_idx * 3 + j;
            auto current_hedge = hedges[hedge_idx];
            int v0_idx = F(face_idx, j);
            int v1_idx = F(face_idx, (j + 1) % 3);
            auto v0 = primal_vertices[v0_idx];
            auto v1 = primal_vertices[v1_idx];
            current_hedge->setStartVertex(v0);
            current_hedge->setEndVertex(v1);
            v0->addIncomingHalfEdge(current_hedge);
            current_hedge->setPrimalFace(face);
            face->addHalfEdge(current_hedge);
            int next_idx = face_idx * 3 + ((j + 1) % 3);
            current_hedge->setNextHalfEdge(hedges[next_idx]);
        }
    }
    std::vector<std::vector<int>> vertex_to_hedges(V.rows());

    for (int i = 0; i < num_hedges; i++) {
        auto hedge = hedges[i];
        int start_idx = hedge->getStartVertex()->index;
        vertex_to_hedges[start_idx].push_back(i);
    }

    for (int i = 0; i < num_hedges; i++) {
        auto hedge = hedges[i];
        if (!hedge->flip.expired()) continue;

        int start_idx = hedge->getEndVertex()->index;
        int end_idx = hedge->getStartVertex()->index;

        bool found_opposite = false;
        for (int hedge_idx : vertex_to_hedges[start_idx]) {
            auto candidate = hedges[hedge_idx];
            if (candidate->getEndVertex()->index == end_idx) {
                hedge->setFlipHalfEdge(candidate);
                candidate->setFlipHalfEdge(hedge);
                found_opposite = true;
                break;
            }
        }

        if (!found_opposite) {
            hedge->boundary = true;
        }
    }
}

//------------------------------------


//void Mesh::vertexDegreeStatistics() {
//
//    std::cout << "Computing vertex degree distribution " << std::endl;
//    auto start = std::chrono::high_resolution_clock::now(); // for measuring time performances
//    int *vDS = new int[V.rows()];
//    //TODO: Implement the vertex degree statistics
//
//    for (int i = 3; i < V.rows(); i++)
//    	if (vDS[i] != 0)
//    		std::cout << "number of degrees = " << i << " ; number of occurences = " << vDS[i] << std::endl;
//    std::cout << "Zero for the other vertex degrees between d=3 and n-1" << std::endl;
//    auto finish = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = finish - start;
//}
void Mesh::vertexDegreeStatistics() {
    std::cout << "Computing vertex degree distribution " << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<int> vertexDegrees(V.rows(), 0);
    for (const auto& face : primal_faces) {
        for (int vertexIndex : face->indices_vertices) {
            vertexDegrees[vertexIndex]++;
        }
    }
    std::cout << "Vertex degree distribution:" << std::endl;
    std::unordered_map<int, int> degreeCount;
    for (int degree : vertexDegrees) {
        degreeCount[degree]++;
    }
    for (const auto& entry : degreeCount) {
        std::cout << "Degree: " << entry.first << ", Count: " << entry.second << std::endl;
    }

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time for vertex degree statistics: " << elapsed.count() << " seconds" << std::endl;
    bool consistent = true;
    for (int i = 0; i < V.rows(); i++) {
        int degreeFromHalfEdges = primal_vertices[i]->incoming_half_edges.size();
        if (vertexDegrees[i] != degreeFromHalfEdges) {
            std::cerr << "Inconsistency found at vertex " << i
                      << ": Face-vertex degree = " << vertexDegrees[i]
                      << ", HalfEdge degree = " << degreeFromHalfEdges << std::endl;
            consistent = false;
        }
    }

    if (consistent) {
        std::cout << "Vertex degrees are consistent with the HalfEdge data structure." << std::endl;
    } else {
        std::cout << "Inconsistencies found in vertex degrees!" << std::endl;
    }
}


Eigen::VectorXd Mesh::compute_gaussian_curvature() {
    Eigen::VectorXd gaussian_curvature(primal_vertices.size());
    gaussian_curvature.setZero();

    for (auto& vertex : primal_vertices) {
        double angle_sum = 0.0;
        double area_sum = 0.0;

        auto one_ring_faces = vertex->getOneRingFaces();
        for (auto& face : one_ring_faces) {
            auto hedge = face->getHalfEdges()[0];
            Eigen::Vector3d v0 = hedge->getStartVertex()->position;
            Eigen::Vector3d v1 = hedge->getNextHalfEdge()->getStartVertex()->position;
            Eigen::Vector3d v2 = hedge->getNextHalfEdge()->getNextHalfEdge()->getStartVertex()->position;

            Eigen::Vector3d e1 = (v1 - v0).normalized();
            Eigen::Vector3d e2 = (v2 - v0).normalized();
            double angle = acos(e1.dot(e2));
            angle_sum += angle;
            Eigen::Vector3d cross = (v1 - v0).cross(v2 - v0);
            double area = cross.norm() / 2.0;
            area_sum += area / 3.0;
        }

        gaussian_curvature(vertex->index) = (2.0 * M_PI - angle_sum) / area_sum;
    }

    return gaussian_curvature;
}


void Mesh::compute_voronoi_area(){
    //TODO
}



void Mesh::count_boundaries() {
    std::unordered_set<int> visited;
    int boundary_count = 0;
    for (const auto& hedge : hedges) {
        if (hedge->boundary && visited.find(hedge->index) == visited.end()) {
            boundary_count++;
            auto current = hedge;

            do {
                visited.insert(current->index);
                current = current->getNextHalfEdge();
            } while (current != hedge);
        }
    }

    std::cout << "boundary_count  " << boundary_count << std::endl;
}



MatrixXd Mesh::compute_vertex_normals(){
    Eigen::MatrixXd vertex_normals(primal_vertices.size(), 3);
    vertex_normals.setZero();

    auto face_normals = compute_face_normals();

    for (size_t i = 0; i < primal_faces.size(); ++i) {
        auto& face = primal_faces[i];
        Eigen::Vector3d face_normal = face_normals.row(i);

        for (auto& vertex : face->getVertices()) {
            vertex_normals.row(vertex->index) += face_normal;
        }
    }

    for (size_t i = 0; i < primal_vertices.size(); ++i) {
        vertex_normals.row(i).normalize();
    }

    return vertex_normals;
}

MatrixXd Mesh::compute_vertex_normals_hed(){
    Eigen::MatrixXd vertex_normals(primal_vertices.size(), 3);
    vertex_normals.setZero();

    for (size_t i = 0; i < primal_vertices.size(); ++i) {
        auto vertex = primal_vertices[i];
        auto one_ring_faces = vertex->getOneRingFaces();

        Eigen::Vector3d norms(0, 0, 0);

        for (auto face : one_ring_faces) {
            auto hedge = face->getHalfEdges()[0];

            Eigen::Vector3d v0 = hedge->getStartVertex()->position;
            Eigen::Vector3d v1 = hedge->getNextHalfEdge()->getStartVertex()->position;
            Eigen::Vector3d v2 = hedge->getNextHalfEdge()->getNextHalfEdge()->getStartVertex()->position;

            Eigen::Vector3d face_normal = (v1 - v0).cross(v2 - v0);
            norms += face_normal.normalized();
        }
        vertex_normals.row(i) = -norms.normalized();
    }

    return vertex_normals;
}

MatrixXd Mesh::compute_face_normals(){

    Eigen::MatrixXd face_normals(primal_faces.size(), 3);
    for (size_t i = 0; i < primal_faces.size(); ++i) {
        auto& face = primal_faces[i];
        auto v0 = face->getVertices()[0]->position;
        auto v1 = face->getVertices()[1]->position;
        auto v2 = face->getVertices()[2]->position;

        Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        n.normalize();
        face_normals.row(i) = n;
    }

    return face_normals;
}

MatrixXd Mesh::compute_face_normals_hed(){
    Eigen::MatrixXd face_normals(primal_faces.size(), 3);
    for (size_t i = 0; i < primal_faces.size(); ++i) {
        auto face = primal_faces[i];
        auto hedge = face->getHalfEdges()[0];

        Eigen::Vector3d v0 = hedge->getStartVertex()->position;
        Eigen::Vector3d v1 = hedge->getNextHalfEdge()->getStartVertex()->position;
        Eigen::Vector3d v2 = hedge->getNextHalfEdge()->getNextHalfEdge()->getStartVertex()->position;

        Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        n.normalize();
        face_normals.row(i) = n;
    }

    return face_normals;
}

