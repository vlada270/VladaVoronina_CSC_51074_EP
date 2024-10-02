#include "TrivConnMesh.h"

#include <iostream>
#include "Eigen/Geometry"
#include "igl/gaussian_curvature.h"
#include "igl/massmatrix.h"
#include <typeinfo>
#include <chrono>


using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::VectorXd;

TrivConnMesh::TrivConnMesh(Eigen::MatrixXd V, Eigen::MatrixXi F) : Mesh(V,F){

    
    
    
    auto start = std::chrono::high_resolution_clock::now();
    initialize_edges();
    auto stop1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = stop1 - start;
    std::cout << "Elapsed time for initializing edges : " << elapsed.count() << " seconds" << std::endl;
    compute_spanning_tree();
    auto stop2 = std::chrono::high_resolution_clock::now();
    elapsed = stop2 - stop1;
    std::cout << "Elapsed time for spanning tree: " << elapsed.count() << " seconds" << std::endl;
    set_frame_field();
    auto stop3 = std::chrono::high_resolution_clock::now();
    elapsed = stop3 - stop2;
    std::cout << "Elapsed time for setting frame field: " << elapsed.count() << " seconds" << std::endl;
    compute_hinge_connection();
    auto stop4 = std::chrono::high_resolution_clock::now();
    elapsed = stop4 - stop3;
    std::cout << "Elapsed time for hinge connection: " << elapsed.count() << " seconds" << std::endl;
    compute_angles_defects();
    auto stop5 = std::chrono::high_resolution_clock::now();
    elapsed = stop5 - stop4;
    std::cout << "Elapsed time for computing angle defects: " << elapsed.count() << " seconds" << std::endl;
    build_LHS();
    auto stop6 = std::chrono::high_resolution_clock::now();
    elapsed = stop6 - stop5;
    std::cout << "Elapsed time for LHS: " << elapsed.count() << " seconds" << std::endl;
    build_RHS_for_problem();
    auto stop7 = std::chrono::high_resolution_clock::now();
    elapsed = stop7 - stop6;
    std::cout << "Elapsed time for RHS: " << elapsed.count() << " seconds" << std::endl;
    solve_for_angles();
    auto stop8 = std::chrono::high_resolution_clock::now();
    elapsed = stop8 - stop7;
    std::cout << "Elapsed time for solve: " << elapsed.count() << " seconds" << std::endl;
    compute_transported_vector_field();
    auto stop9 = std::chrono::high_resolution_clock::now();
    elapsed = stop9 - stop8;
    std::cout << "Elapsed time for parallel transport: " << elapsed.count() << " seconds" << std::endl;

};

TrivConnMesh::TrivConnMesh() : Mesh(){};

double TrivConnMesh::adjustToRange(double angle) {
    const double twoPi = 2 * M_PI;
    const double pi = M_PI;
    // Ensure the angle is within the range [0, 2*pi]
    while (angle < -pi) {
        angle += twoPi;
    }
    while (angle > pi) {
        angle -= twoPi;
    }
    return angle;
};

void TrivConnMesh::set_frame_field(){

    //TODO
}

void TrivConnMesh::initialize_edges(){
    //TODO
}

/**
 * Computes the hinge connection angles for each edge in the mesh.
 * The hinge connection angle is the angle between the frames of the adjacent faces
 * connected by the edge.
 * 
 * This function iterates over each edge in the mesh and calculates the hinge connection
 * angle using the frames of the adjacent faces. The angle is then assigned to the edge
 * and its corresponding half-edge.
 */
void TrivConnMesh::compute_hinge_connection()
{
    //TODO
}

/**
 * Computes the angles and defects for each vertex in the mesh.
 * 
 * This method iterates over the primal vertices of the mesh and computes the angles and defects for each vertex.
 * The angles and defects are important measures in the field of computational geometry and finite element analysis.
 * 
 * @note This method assumes that the primal vertices of the mesh have already been initialized.
 * 
 * @see Mesh::primal_vertices
 */
void TrivConnMesh::compute_angles_defects(){

  //TODO
}

void TrivConnMesh::compute_spanning_tree(){

  //TODO

}

void TrivConnMesh::build_LHS(){
    //TODO
    

};

void TrivConnMesh::build_RHS_for_problem(){
    singularities = Eigen::VectorXd::Zero(V.rows());
    double pi = M_PI;
    // singularities(120) = 2.;

    //TODO

    check_RHS_for_consistency(singularities);

};

void TrivConnMesh::check_RHS_for_consistency(const VectorXd& singularities){

    int sum_singularities = 0;
    for(int i = 0; i<singularities.size(); i++){
        sum_singularities+=singularities(i);
    }
    // std::cout<<" Euler Characteristic "<<this->euler_characteristic<<" sum singularities "<< sum_singularities<<"\n";


    if(sum_singularities==this->euler_characteristic){
        std::cout<<"Go ahead, these singularities make sense!"<<std::endl;
    }else{
        std::cout<<" Wrong number of singularities \n";
        std::cout<<" Euler Characteristic "<<this->euler_characteristic<<" sum singularities "<< sum_singularities<<"\n";
        std::cout<<" The vector field will have singularities!!!"<<std::endl;
    }
    
};

// solve for the adjustment angles.
void TrivConnMesh::solve_for_angles(){
    //TODO

};

//now propagate the vector to obtain the field
void TrivConnMesh::compute_transported_vector_field(){

    //TODO

};


