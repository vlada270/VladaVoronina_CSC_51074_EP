// Include necessary headers
#include "gtest/gtest.h"
#include <iostream>
#include <Eigen/Core>
#include "igl/readPLY.h"
#include "../src/cubic.cpp"
#include "../src/lagrange.cpp"
#include "../src/linear.cpp"
#include "../src/hermite.cpp"
#include "../src/util.h"

Eigen::MatrixXd V1;
Eigen::MatrixXi F1;

TEST(Lagrange_check,check_interpolation){

  #ifdef USING_GNU
            // std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
            igl::readPLY("../data/curve.ply", V1, F1);
            std::string filename = "../data/curve.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load GNU-specific meshes
    #elif defined(USING_MSVC)
            // std::cout << "Using MSVC compiler. Loading MSVC-specific meshes." << std::endl;
            igl::readPLY("../../../../data/curve.ply", V1, F1);
            std::string filename = "../../../../../data/curve.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load MSVC-specific meshes
    #else
            // std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
            igl::readPLY("../data/curve.ply", V, F);
            std::string filename = "../data/curve.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load default meshes
    #endif  
    LagrangeInterpolation interp_lagrange = LagrangeInterpolation(V1);


    int resolution = 37;
    Eigen::MatrixXd linspace = MatrixXd::Zero(resolution, 3);
    util::build_linspace(linspace, V1); // initialize the X axis of the interpolation
    std::vector<Eigen::Vector3d> points;
    std::vector<std::array<int, 2>> edges;
    std::vector<std::array<double, 3>> colors;

    for (size_t i = 0; i < resolution; i++) {
    double time = linspace(i, 0);
    linspace(i, 1) = interp_lagrange.eval_function(time);
    Eigen::Vector3d pt = Eigen::Vector3d(linspace(i, 0), linspace(i, 1), 0);
    points.push_back(pt);
    if(i<resolution-1){
        edges.push_back({int(i), int(i+1)});
        colors.push_back({0.0, 0.0, 1.0});
    }
    }

    Eigen::MatrixXd row_after(1,3);
    
    row_after(0,0) = -0.686667;
    row_after(0,1) = -0.343917;
    row_after(0,2) = 0.0;

    ASSERT_NEAR(row_after(0,0),linspace(12,0),0.0001);
    ASSERT_NEAR(row_after(0,1),linspace(12,1),0.0001);
    ASSERT_NEAR(row_after(0,2),linspace(12,2),0.0001);
}


TEST(Cubic_check,check_interpolation){

  #ifdef USING_GNU
            // std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
            igl::readPLY("../data/curve1.ply", V1, F1);
            std::string filename = "../data/curve1.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load GNU-specific meshes
    #elif defined(USING_MSVC)
            // std::cout << "Using MSVC compiler. Loading MSVC-specific meshes." << std::endl;
            igl::readPLY("../../../../data/curve1.ply", V1, F1);
            std::string filename = "../../../../../data/curve1.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load MSVC-specific meshes
    #else
            // std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
            igl::readPLY("../data/curve1.ply", V, F);
            std::string filename = "../data/curve.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load default meshes
    #endif  
    CubicInterpolation interp_cubic = CubicInterpolation(V1);


    int resolution = 49;
    Eigen::MatrixXd linspace = MatrixXd::Zero(resolution, 3);
    util::build_linspace(linspace, V1); // initialize the X axis of the interpolation
    std::vector<Eigen::Vector3d> points;
    std::vector<std::array<int, 2>> edges;
    std::vector<std::array<double, 3>> colors;

    for (size_t i = 0; i < resolution; i++) {
    double time = linspace(i, 0);
    linspace(i, 1) = interp_cubic.eval_function(time);
    Eigen::Vector3d pt = Eigen::Vector3d(linspace(i, 0), linspace(i, 1), 0);
    points.push_back(pt);
    if(i<resolution-1){
        edges.push_back({int(i), int(i+1)});
        colors.push_back({0.0, 0.0, 1.0});
    }
    }

    Eigen::MatrixXd row_after(1,3);
    row_after(0,0) = -0.40625;
    row_after(0,1) = 1.15878;
    row_after(0,2) = 0.0;
    Eigen::MatrixXd linspace_row = linspace.row(17);
    std::cout<<"linspace cubic"<<linspace_row<<std::endl;
    

    ASSERT_NEAR(row_after(0,0),linspace_row(0,0),0.0001);
    ASSERT_NEAR(row_after(0,1),linspace_row(0,1),0.0001);
    ASSERT_NEAR(row_after(0,2),linspace_row(0,2),0.0001);
}



TEST(Linear_check,check_interpolation){

  #ifdef USING_GNU
            // std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
            igl::readPLY("../data/curve.ply", V1, F1);
            // std::string filename = "../data/curve1.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load GNU-specific meshes
    #elif defined(USING_MSVC)
            // std::cout << "Using MSVC compiler. Loading MSVC-specific meshes." << std::endl;
            igl::readPLY("../../../../data/curve.ply", V1, F1);
            // std::string filename = "../../../../../data/curve1.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load MSVC-specific meshes
    #else
            // std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
            igl::readPLY("../data/curve.ply", V, F);
            // std::string filename = "../data/curve.ply";
            // std::cout << "loading: " << filename << std::endl;
            // Load default meshes
    #endif  
    LinearInterpolation interp_linear = LinearInterpolation(V1);


    int resolution = 69;
    Eigen::MatrixXd linspace = MatrixXd::Zero(resolution, 3);
    util::build_linspace(linspace, V1); // initialize the X axis of the interpolation
    std::vector<Eigen::Vector3d> points;
    std::vector<std::array<int, 2>> edges;
    std::vector<std::array<double, 3>> colors;

    for (size_t i = 0; i < resolution; i++) {
    double time = linspace(i, 0);
    linspace(i, 1) = interp_linear.eval_function(time);
    Eigen::Vector3d pt = Eigen::Vector3d(linspace(i, 0), linspace(i, 1), 0);
    points.push_back(pt);
    if(i<resolution-1){
        edges.push_back({int(i), int(i+1)});
        colors.push_back({0.0, 0.0, 1.0});
    }
    }

    Eigen::MatrixXd row_after(1,3);
    row_after(0,0) = 0.13;
    row_after(0,1) = 0.316154;
    row_after(0,2) = 0.0;

    ASSERT_NEAR(row_after(0,0),linspace(34,0),0.0001);
    ASSERT_NEAR(row_after(0,1),linspace(34,1),0.0001);
    ASSERT_NEAR(row_after(0,2),linspace(34,2),0.0001);
}