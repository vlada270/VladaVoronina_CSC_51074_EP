// Include necessary headers
#include "gtest/gtest.h"
#include <iostream>
#include <Eigen/Core>
#include "igl/readOBJ.h"
#include "igl/readOFF.h"
#include "../src/Transform.h"


Eigen::MatrixXd V0;
Eigen::MatrixXi F0;
using Eigen::MatrixXd, Eigen::Vector3d, Eigen::RowVector3d, Eigen::Quaterniond, Eigen::MatrixXi;



TEST(TransformCheck,BasicSanityCheck){
    #ifdef USING_GNU
        int i0 = igl::readOFF("../data/star.off", V0, F0);
        // Load GNU-specific meshes
    #elif defined(USING_MSVC)
        int i0 = igl::readOFF("../../../../../data/star.off", V0, F0);
        // Load MSVC-specific meshes
    #else
        int i0 = igl::readOFF("../data/star.off", V0, F0);
        // Load default meshes
    #endif

    MatrixXd first_vertex = V0.row(0);
    first_vertex = 3*first_vertex;

    Transform T;

    T.uniform_scale(3.0);

    T.apply_transform(V0);

    MatrixXd transformed_vertex = V0.row(0);

    ASSERT_EQ(first_vertex,transformed_vertex);
}

TEST(TransformCheck,Check_Eq_operator){
    #ifdef USING_GNU
        int i0 = igl::readOFF("../data/homer.off", V0, F0);
        // Load GNU-specific meshes
    #elif defined(USING_MSVC)
        int i0 = igl::readOFF("../../../../../data/homer.off", V0, F0);
        // Load MSVC-specific meshes
    #else
        int i0 = igl::readOFF("../data/homer.off", V0, F0);
        // Load default meshes
    #endif

    Transform T1;
    Transform T2;

    T1.uniform_scale(3.0);
    T2.uniform_scale(2.0);
    ASSERT_EQ((T1*T2).M,(T2*T1).M);
    ASSERT_EQ((T1*T2).M,(T1.compose(T2)).M);
    ASSERT_EQ(T1.M*T2.M,(T1*T2.M).M);
    ASSERT_EQ(T1.M*T2.M,(T1.compose(T2.M)).M);
}

TEST(TransformCheck, Transform_Composition_Check){
    #ifdef USING_GNU
        int i0 = igl::readOFF("../data/star.off", V0, F0);
        // Load GNU-specific meshes
    #elif defined(USING_MSVC)
        int i0 = igl::readOFF("../../../../../data/star.off", V0, F0);
        // Load MSVC-specific meshes
    #else
        int i0 = igl::readOFF("../data/star.off", V0, F0);
        // Load default meshes
    #endif
	Transform T;
	T.uniform_scale(3.0);
	T.apply_transform(V0);
	// std::cout<<"V after transformation: \n"<<V<<std::endl;
	Transform T2;
	T2.translate(Eigen::Vector3d(1.7, 1.0, 14.6));
	T2.apply_transform(V0);
	

    MatrixXd row_after(1,3);
    row_after(0,0) = 1.7;
    row_after(0,1) =  13.0;
    row_after(0,2) =  14.6;
    
    
    MatrixXd transformed_row = V0.row(9);
	ASSERT_NEAR(row_after(0),transformed_row(0),0.0001);
    ASSERT_NEAR(row_after(1),transformed_row(1),0.0001);
    ASSERT_NEAR(row_after(2),transformed_row(2),0.0001);
}


TEST(TransformCheck, QuaternionCheck){
    #ifdef USING_GNU
        int i0 = igl::readOFF("../data/homer.off", V0, F0);
        // Load GNU-specific meshes
    #elif defined(USING_MSVC)
        int i0 = igl::readOFF("../../../../../data/homer.off", V0, F0);
        // Load MSVC-specific meshes
    #else
        int i0 = igl::readOFF("../data/homer.off", V0, F0);
        // Load default meshes
    #endif
    Transform T;
	T.rotate_with_quaternions(V0, Eigen::Vector3d(1.0, 2.1, 0.1), 3.14/3);
    MatrixXd row_after(1,3);
    row_after(0,0) = 0.490936;
    row_after(0,1) =  0.185639;
    row_after(0,2) =  -0.0898287;

    MatrixXd transformed_row = V0.row(9);
    ASSERT_NEAR(row_after(0),transformed_row(0),0.0001);
    ASSERT_NEAR(row_after(1),transformed_row(1),0.0001);
    ASSERT_NEAR(row_after(2),transformed_row(2),0.0001);

}