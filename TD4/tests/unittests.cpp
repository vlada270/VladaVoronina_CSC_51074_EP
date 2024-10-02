#include "gtest/gtest.h"
#include "Eigen/Core"

#include "../src/Mesh.h"
#include "igl/readOBJ.h"
#include "igl/readOFF.h"
#include "igl/gaussian_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"

Eigen::MatrixXd V0;
Eigen::MatrixXi F0;

#ifdef USING_GNU
    int i0 = igl::readOFF("../data/bunny_new.off", V0, F0);
    // Load GNU-specific meshes
#elif defined(USING_MSVC)
    int i0 = igl::readOFF("../../../../../data/bunny_new.off", V0, F0);
    // Load MSVC-specific meshes
#else
    int i0 = igl::readOFF("../data/bunny_new.off", V0, F0);
    // Load default meshes
#endif


Mesh test_mesh = Mesh(V0,F0);


TEST(MeshTest, CheckBasicSanityMesh)
{
    EXPECT_EQ(test_mesh.primal_vertices.size(),V0.rows() );
    EXPECT_EQ(test_mesh.primal_faces.size(),F0.rows() );
}


// Test that the pointer is initially null
TEST(VertexTest, CheckGoodInitialization) {

    int n_vertices = test_mesh.primal_vertices.size();

    for(int idx_vtx = 0; idx_vtx<n_vertices; idx_vtx++){
        
        VertexPtr vtx = test_mesh.primal_vertices[idx_vtx];
        std::vector<PrimalFacePtr> adjacent_faces = vtx->getOneRingFaces();
        std::vector<HalfEdgePtr> incoming_hedges = vtx->getIncomingHalfEdges();
        EXPECT_GT(adjacent_faces.size(), 0);
        EXPECT_GT(incoming_hedges.size(), 0);
        for (PrimalFacePtr ptr : adjacent_faces) {
            ASSERT_NE(ptr, nullptr);
        }

        for (HalfEdgePtr ptr : incoming_hedges) {
            ASSERT_NE(ptr, nullptr);
            ASSERT_NE(ptr->getStartVertex(), nullptr);
            ASSERT_NE(ptr->getEndVertex(), nullptr);
            ASSERT_NE(ptr->getFlipHalfEdge(), nullptr);
            ASSERT_NE(ptr->getNextHalfEdge(), nullptr);
            ASSERT_NE(ptr->getPrimalFace(), nullptr);
        }

        EXPECT_LE(vtx->index,n_vertices);

    
    }


    
}

// Test case for verifying if all the pointer attributes are indeed attributed
TEST(HalfEdgeTest, CheckGoodInitialization) {
    
    for(int hedge_itr = 0; hedge_itr<test_mesh.hedges.size(); hedge_itr++)
    {   
        HalfEdgePtr hedge = test_mesh.hedges[hedge_itr];
        //make sure that the pointer attributes are actually all initialized
        ASSERT_NE(hedge->getStartVertex(), nullptr);
        ASSERT_NE(hedge->getEndVertex(), nullptr);
        ASSERT_NE(hedge->getFlipHalfEdge(), nullptr);
        ASSERT_NE(hedge->getNextHalfEdge(), nullptr);
        ASSERT_NE(hedge->getPrimalFace(), nullptr);

        ASSERT_NE(hedge->getFlipHalfEdge()->getStartVertex(), hedge->getEndVertex());
        ASSERT_NE(hedge->getFlipHalfEdge()->getEndVertex(), hedge->getStartVertex());
        ASSERT_NE(hedge->getFlipHalfEdge()->getFlipHalfEdge(), hedge);
        ASSERT_NE(hedge->getFlipHalfEdge()->getNextHalfEdge(), nullptr);
        ASSERT_NE(hedge->getFlipHalfEdge()->getPrimalFace(), nullptr);

    }

}


// Test case for verifying if the two vectors contain the same pointers for each Vertex
TEST(HalfEdgeTest, SanityCheckFlipFlip) {
        
    //Check if the flip attribute is set properly
    for(int hedge_itr = 0; hedge_itr<test_mesh.hedges.size(); hedge_itr++)
    {
        HalfEdgePtr hedge = test_mesh.hedges[hedge_itr];
        EXPECT_EQ(hedge, (hedge->getFlipHalfEdge()->getFlipHalfEdge()));
    }

}

//Test to make sure that everywhere there is a pointer as a class attribute, it is indeed not a nullpointer
TEST(HalfEdgeTest, SanityCheckNext) {
    
    //Check if for a triangle mesh there is actually next next next fine
    
    for(int hedge_itr = 0; hedge_itr<test_mesh.hedges.size(); hedge_itr++)
    {
        HalfEdgePtr hedge = test_mesh.hedges[hedge_itr];
        PrimalFacePtr face = hedge->getPrimalFace();
        EXPECT_EQ(hedge, (hedge->getNextHalfEdge()->getNextHalfEdge()->getNextHalfEdge()));
        EXPECT_EQ(face, (hedge->getNextHalfEdge()->getPrimalFace()));
        EXPECT_EQ(face, (hedge->getNextHalfEdge()->getNextHalfEdge()->getPrimalFace()));
        EXPECT_EQ(face, (hedge->getNextHalfEdge()->getNextHalfEdge()->getNextHalfEdge()->getPrimalFace()));
    }

}

TEST(HalfEdgeTest, SanityCheckPrimalFace) {
    
    //Check if the pointer attributes are actually all initialized
    for(int hedge_itr = 0; hedge_itr<test_mesh.hedges.size(); hedge_itr++)
    {
        HalfEdgePtr hedge = test_mesh.hedges[hedge_itr];
        ASSERT_NE(hedge->getPrimalFace(), nullptr);
        ASSERT_NE(hedge->getFlipHalfEdge()->getPrimalFace(), nullptr);
        ASSERT_NE(hedge->getNextHalfEdge()->getPrimalFace(), nullptr);
        ASSERT_EQ(hedge->getPrimalFace(),hedge->getFlipHalfEdge()->getFlipHalfEdge()->getPrimalFace());
    }

}


TEST(LaplacianTest,CheckWithLibigl) {
    //TODO
    // Use the cotmatrix from libigl to compare with the dirichlet matrix from the mesh
}

TEST(AreaMatrixTest, CheckWithLibigl) {
    //TODO
    // Use the massmatrix from libigl to make sure that the calculation of the voronoi area works
}
