// Include necessary headers
#include "gtest/gtest.h"
#include <iostream>
#include <Eigen/Core>
#include "igl/readPLY.h"
#include "../src/Bezier.h"
#include "../src/Util.h"

Eigen::MatrixXd V1;
Eigen::MatrixXi F1;

TEST(DeCasteljeauCheck,check_de_casteljeau){

  #ifdef USING_GNU
            // std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
            igl::readPLY("../data/hair.ply", V1, F1);
            std::string filename = "../data/hair.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load GNU-specific meshes
    #elif defined(USING_MSVC)
            // std::cout << "Using MSVC compiler. Loading MSVC-specific meshes." << std::endl;
            igl::readPLY("../../../../data/hair.ply", V1, F1);
            std::string filename = "../../../../../data/hair.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load MSVC-specific meshes
    #else
            // std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
            igl::readPLY("../data/hair.ply", V1, F1);
            std::string filename = "../data/hair.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load default meshes
    #endif  
    
    Bezier bezier;
    int resolution = 162;
    Eigen::MatrixXd test_value = bezier.de_casteljau(V1,0.6432);
    double test_value_1 = test_value(0,0);
    double test_value_2 = test_value(0,1);
    double test_value_3 = test_value(0,2);
    Eigen::MatrixXd row_after(1,3);
    row_after(0,0) = 5.44681;
    row_after(0,1) = -12.0808;
    row_after(0,2) = 1.1869;
    ASSERT_NEAR(row_after(0,0),test_value(0,0),0.0001);
    ASSERT_NEAR(row_after(0,1),test_value(0,1),0.0001);
    ASSERT_NEAR(row_after(0,2),test_value(0,2),0.0001);

}