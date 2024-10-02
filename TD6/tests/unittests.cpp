#include "gtest/gtest.h"
#include "Eigen/Core"

#include "../src/Conformal.h"
#include "igl/readOBJ.h"
#include "igl/readOFF.h"
#include "igl/gaussian_curvature.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"


TEST(SanityCheckDirichlet,size_check){
    #ifdef USING_GNU
        int i0 = igl::readOFF("../data/beetle.off", V0, F0);
        // Load GNU-specific meshes
    #elif defined(USING_MSVC)
        int i0 = igl::readOFF("../../../../../data/beetle.off", V0, F0);
        // Load MSVC-specific meshes
    #else
        int i0 = igl::readOFF("../data/beetle.off", V0, F0);
        // Load default meshes
    #endif

    int n_vertices = V0.rows();
    ConformalParametrization parameterization = ConformalParametrization(V0,F0);

    parameterization.build_parametrizations();
    ASSERT_EQ(parameterization.Dirichlet.rows(),2*n_vertices);
    ASSERT_EQ(parameterization.Dirichlet.cols(),2*n_vertices);
    ASSERT_EQ(parameterization.Area.rows(),2*n_vertices);
    ASSERT_EQ(parameterization.Area.cols(),2*n_vertices);
    ASSERT_EQ(parameterization.ConformalEnergy.rows(),2*n_vertices);
    ASSERT_EQ(parameterization.ConformalEnergy.cols(),2*n_vertices);

}


TEST(AreaCheck, AreaCheck_Ar_Test){
    //TODO
    // Build yourself a 2D parameterization of a square with a sidelength of 1.68
    // Put for instance one vertex at (0.0) and then use 100 triangles to build a triangulation of the square
    // Then check that the area, that you can compute by hand matches 
    // flattened_vector_square.transpose()*Area*flattened_vector_square

}

