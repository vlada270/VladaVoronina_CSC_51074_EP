/**
 * @file unittest.cpp
 * @brief Unit tests for the Mesh class and its related classes and functions.
 *
 * This file contains unit tests for the Mesh class and its related classes and functions.
 * The tests cover various aspects of the Mesh class, including initialization, pointer attributes,
 * and consistency checks. The tests verify if the Mesh class behaves as expected and if the
 * pointer attributes are properly initialized and connected.
 *
 * The unit tests are implemented using the Google Test framework and cover different aspects
 * of the Mesh class, such as vertex initialization, half-edge initialization, edge initialization,
 * primal face initialization, and connection tests. The tests check if the pointer attributes
 * of the Mesh class and its related classes are not null and if they are properly connected.
 * Additionally, the tests verify if the angle defect of the hinge connection corresponds to
 * the Gaussian curvature of the mesh.
 *
 * The unit tests are organized into different test cases, each focusing on a specific aspect
 * of the Mesh class. The test cases include:
 * - MeshTest: Checks the basic sanity of the Mesh class, including the size of primal vertices
 *   and primal faces.
 * - VertexTest: Verifies the good initialization of vertex pointer attributes and checks if
 *   the vertex index is within the expected range.
 * - HalfEdgeTest: Tests the good initialization of half-edge pointer attributes, checks the
 *   sanity of half-edges with edges, and verifies the consistency of the flip attribute.
 * - EdgeTest: Checks the good initialization of edge pointer attributes and verifies the
 *   consistency of edges with half-edges.
 * - PrimalFaceTest: Verifies the good initialization of primal face pointer attributes and
 *   checks the presence of edges, half-edges, vertices, and indices.
 * - ConnectionTest: Tests the connection between the angle defect of the hinge connection
 *   and the Gaussian curvature of the mesh.
 *
 * The unit tests are implemented using assertions from the Google Test framework. The tests
 * use various assertions, such as EXPECT_EQ, EXPECT_GT, ASSERT_NE, EXPECT_NEAR, to check
 * the expected behavior and verify the correctness of the Mesh class and its related classes
 * and functions.
 */
// Include necessary headers
#include "gtest/gtest.h"
#include "Eigen/Core"
// #include "../src/Mesh.h" // Assuming this header includes the definition of Mesh
#include "../src/TrivConnMesh.h"

#include "igl/readOBJ.h"
#include "igl/readOFF.h"
#include "igl/gaussian_curvature.h"


Eigen::MatrixXd V0;
Eigen::MatrixXi F0;

#ifdef USING_GNU
    int i0 = igl::readOBJ("../data/spot.obj", V0, F0);
    // Load GNU-specific meshes
#elif defined(USING_MSVC)
    int i0 = igl::readOBJ("../../../../../data/spot.obj", V0, F0);
    // Load MSVC-specific meshes
#else
    int i0 = igl::readOBJ("../data/spot.obj", V0, F0);
    // Load default meshes
#endif
//int i0 = igl::readOBJ("../data/spot.obj",V0,F0);
// int i0 = igl::readOBJ("../../../../../data/spot.obj", V0, F0);

// int i0 = igl::readOFF("../data/torus_fine.off",V0,F0);

// Mesh test_mesh = Mesh(V0,F0);
TrivConnMesh test_mesh = TrivConnMesh(V0,F0);


TEST(MeshTest, CheckBasicSanityMesh)
{
    EXPECT_EQ(test_mesh.primal_vertices.size(),V0.rows() );
    EXPECT_EQ(test_mesh.primal_faces.size(),F0.rows() );
    EXPECT_EQ(2*test_mesh.edges.size(),test_mesh.hedges.size());
}

TEST(HedgeTest, CheckInit){
    for(int i = 0; i<test_mesh.hedges.size(); i++){
        HalfEdgePtr current_hedge = test_mesh.hedges[i];
        ASSERT_NE(current_hedge, nullptr);
        ASSERT_NE(current_hedge->getEdge(), nullptr);
        ASSERT_NE(current_hedge->getEndVertex(), nullptr);
        ASSERT_NE(current_hedge->getStartVertex(), nullptr);
        ASSERT_NE(current_hedge->getPrimalFace(), nullptr);
    }

}

TEST(EdgeTest, CheckInit){
    for(int i = 0; i<test_mesh.edges.size(); i++){
        EdgePtr current_edge = test_mesh.edges[i];
        ASSERT_NE(current_edge, nullptr);
        ASSERT_NE(current_edge->getHedge(), nullptr);
        ASSERT_NE(current_edge->getLeftFace(), nullptr);
        ASSERT_NE(current_edge->getRightFace(), nullptr);
        ASSERT_NE(current_edge->getStartVertex(), nullptr);
        ASSERT_NE(current_edge->getEndVertex(), nullptr);
    }
}

TEST(ConnectionTest, CheckGaussianCurvature) {
    double tol = 1e-6;

    Eigen::MatrixXd K = Eigen::MatrixXd::Ones(test_mesh.V.rows(),1);
    
    igl::gaussian_curvature(test_mesh.V,test_mesh.F,K);
    // for(auto vertex_itr = test_mesh.primal_vertices.begin(); vertex_itr!=test_mesh.primal_vertices.end(); vertex_itr++){
    for(int i = 0; i<test_mesh.primal_vertices.size(); i++){
        VertexPtr vertex_itr = test_mesh.primal_vertices[i];
        double angle_defect = vertex_itr->angle_defect;
        double gaussian_curvature = K(vertex_itr->index);
        int index = round((double(angle_defect) - double(K(vertex_itr->index,0)))/(2*M_PI));
        double cleared_defect = angle_defect - 2*M_PI*index;
        EXPECT_NEAR(cleared_defect,gaussian_curvature,tol);
        
    }
}

// Test that the pointer is initially null
TEST(VertexTest, CheckGoodInitialization) {

    int n_vertices = V0.rows();

    for(auto vertex_itr = test_mesh.primal_vertices.begin(); vertex_itr!=test_mesh.primal_vertices.end(); vertex_itr++){

        std::vector<PrimalFace*> adjacent_faces = vertex_itr->one_ring_faces;
        std::vector<HalfEdge*> outgoing_hedges = vertex_itr->outgoing_half_edges;
        EXPECT_GT(adjacent_faces.size(), 0);
        EXPECT_GT(outgoing_hedges.size(), 0);
        for (PrimalFace* ptr : adjacent_faces) {
            ASSERT_NE(ptr, nullptr);
        }

        for (HalfEdge* ptr : outgoing_hedges) {
            ASSERT_NE(ptr, nullptr);
        }

        EXPECT_LE(vertex_itr->index,n_vertices);
    
    }
    
}

// Test case for verifying if all the pointer attributes are indeed attributed
TEST(HalfEdgeTest, CheckGoodInitialization) {
    
    for(auto hedge_itr = test_mesh.hedges.begin(); hedge_itr!=test_mesh.hedges.end(); hedge_itr++)
    {
        //make sure that the pointer attributes are actually all initialized
        ASSERT_NE(hedge_itr->start, nullptr);
        ASSERT_NE(hedge_itr->end, nullptr);
        ASSERT_NE(hedge_itr->flip, nullptr);
        ASSERT_NE(hedge_itr->next, nullptr);
        ASSERT_NE(hedge_itr->primal_face, nullptr);
        ASSERT_NE(hedge_itr->edge, nullptr);
    }

}

// Test case for verifying if all the pointer attributes are indeed attributed
TEST(HalfEdgeTest, CheckSanityWithEdges) {
    int n_hedges = test_mesh.hedges.size();
    int n_edges = test_mesh.edges.size();
    EXPECT_EQ(2*n_edges,n_hedges);
    for(auto hedge_itr = test_mesh.hedges.begin(); hedge_itr!=test_mesh.hedges.end(); hedge_itr++)
    {
        EXPECT_LE(hedge_itr->index,test_mesh.hedges.size());
        if(hedge_itr->sign_edge==1){
            ASSERT_NE(hedge_itr->edge->v_primal_start, nullptr);
            ASSERT_NE(hedge_itr->edge->v_primal_start, nullptr);
            EXPECT_EQ(hedge_itr->start->index,hedge_itr->edge->v_primal_start->index);
            EXPECT_EQ(hedge_itr->end->index,hedge_itr->edge->v_primal_end->index);
        }
        if(hedge_itr->sign_edge==-1){
            ASSERT_NE(hedge_itr->edge->v_primal_start, nullptr);
            ASSERT_NE(hedge_itr->edge->v_primal_start, nullptr);
            EXPECT_EQ(hedge_itr->end->index,hedge_itr->edge->v_primal_start->index);
            EXPECT_EQ(hedge_itr->start->index,hedge_itr->edge->v_primal_end->index);   
        }
    }
    

}

// Test case for verifying if the two vectors contain the same pointers for each Vertex
TEST(HalfEdgeTest, SanityCheckFlipFlip) {
        
    //Check if the flip attribute is set properly
    for( auto hedge_itr = test_mesh.hedges.begin(); hedge_itr!=test_mesh.hedges.end(); hedge_itr++)
    {
        EXPECT_EQ(*hedge_itr, *(hedge_itr->flip->flip));
    }

}

//Test to make sure that everywhere there is a pointer as a class attribute, it is indeed not a nullpointer
TEST(HalfEdgeTest, SanityCheckNext) {
    
    //Check if for a triangle mesh there is actually next next next fine
    
    for( auto hedge_itr = test_mesh.hedges.begin(); hedge_itr!=test_mesh.hedges.end(); hedge_itr++)
    {
        EXPECT_EQ(*hedge_itr, *(hedge_itr->next->next->next));
        EXPECT_EQ(*(hedge_itr->edge), *(hedge_itr->flip->edge));
    }

}

TEST(EdgeTest, CheckGoodInitialization){
    
    //make sure that the pointer attributes are actually all initialized
    for(auto edge_itr = test_mesh.edges.begin(); edge_itr!=test_mesh.edges.end(); edge_itr++)
    {
        ASSERT_NE(edge_itr->v_primal_start, nullptr);
        ASSERT_NE(edge_itr->v_primal_end, nullptr);
        ASSERT_NE(edge_itr->face_left_ptr, nullptr);
        ASSERT_NE(edge_itr->face_right_ptr, nullptr);
        ASSERT_NE(edge_itr->hedge_with_frame, nullptr);
        ASSERT_NE(edge_itr->hedge, nullptr);
    }

}

TEST(EdgeTest, CheckConsistencyWithHedge){
    
    //make sure that the pointer attributes are actually all initialized
    for(auto edge_itr = test_mesh.edges.begin(); edge_itr!=test_mesh.edges.end(); edge_itr++)
    {
        EXPECT_LE(edge_itr->index,test_mesh.edges.size());
        EXPECT_LE(edge_itr->hedge->index, test_mesh.hedges.size());
        ASSERT_NE(edge_itr->hedge->edge,nullptr);
        EXPECT_EQ(edge_itr->hedge->start->index,edge_itr->v_primal_start->index);
        EXPECT_EQ(edge_itr->hedge->end->index,edge_itr->v_primal_end->index);
        EXPECT_EQ(edge_itr->hedge->flip->start->index,edge_itr->v_primal_end->index);
        EXPECT_EQ(edge_itr->hedge->flip->end->index,edge_itr->v_primal_start->index);
    }

}

TEST(PrimalFaceTest, CheckGoodInitialization) {

    for(auto face_itr = test_mesh.primal_faces.begin(); face_itr!=test_mesh.primal_faces.end(); face_itr++){

        std::vector<Edge*> edges_face = face_itr->edges_face;
        std::vector<HalfEdge*> hedges = face_itr->hedges_face;
        std::vector<Vertex*> vertices_face = face_itr->vertices_face;
        std::vector<int> indices_vertices = face_itr->indices_vertices;
        std::vector<int> indices_edges = face_itr->indices_edges;
        std::vector<int> signs_edges = face_itr->signs_edges;
        
        EXPECT_GT(edges_face.size(), 0);
        EXPECT_GT(hedges.size(), 0);
        EXPECT_GT(vertices_face.size(), 0);
        EXPECT_GT(indices_vertices.size(), 0);
        EXPECT_GT(indices_edges.size(), 0);
        EXPECT_GT(indices_edges.size(), 0);

        EXPECT_EQ(indices_edges.size(),signs_edges.size());
        EXPECT_EQ(indices_edges.size(),edges_face.size());
        EXPECT_EQ(indices_vertices.size(),vertices_face.size());
        
        for (Edge* ptr : edges_face) {
            ASSERT_NE(ptr, nullptr);
        }

        for (HalfEdge* ptr : hedges) {
            ASSERT_NE(ptr, nullptr);
        }

        for (Vertex* ptr : vertices_face) {
            ASSERT_NE(ptr, nullptr);
        }
    }
    
}

//Verify if the angle defect of the hinge connection indeed amounts to the gaussian curvature
TEST(ConnectionTest, CheckGaussianCurvature) {
    double tol = 1e-6;

    Eigen::MatrixXd K = Eigen::MatrixXd::Ones(test_mesh.V.rows(),1);
    
    igl::gaussian_curvature(test_mesh.V,test_mesh.F,K);
    //TODO: Check if the angle defect of the hinge connection indeed amounts to the gaussian curvature up to 2\pi
}





