#include "LaplacianMesh.h"
#include "Eigen/IterativeLinearSolvers"

using Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd, Eigen::Vector3d;


LaplacianMesh::LaplacianMesh(/* args */){}

LaplacianMesh::~LaplacianMesh(){}

void LaplacianMesh::compute_dirichlet(){
    //TODO
}

void LaplacianMesh::compute_area_matrix(){
    //TODO
    
}

void LaplacianMesh::compute_laplacian(){
    //TODO

}


LaplacianMesh::LaplacianMesh(Eigen::MatrixXd V,Eigen::MatrixXi F):Mesh(V,F){
    this->V = V;
    this->F = F;

    compute_dirichlet();

    compute_area_matrix();

    compute_laplacian();

}

void LaplacianMesh::setHeat(Eigen::VectorXd & u){
    //TODO
}

void LaplacianMesh::heat_step_explicit(Eigen::VectorXd & u, double time_step ){
    //TODO
}

void LaplacianMesh::heat_step_implicit(Eigen::VectorXd & u, double time_step){
    //TODO
}