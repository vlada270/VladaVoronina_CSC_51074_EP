#include "ElectricMesh.h"

#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/per_vertex_normals.h"

#include "Eigen/Sparse"

using namespace Eigen;

ElectricMesh::ElectricMesh(MatrixXd V, MatrixXi F):Mesh(V,F){

    set_frame_field();

};

void ElectricMesh::set_frame_field(){

    //TODO
}

void ElectricMesh::initialize_charge_density(const Eigen::MatrixXd rho_in)
{

    //TODO
};

void ElectricMesh::solve_for_u()
{

    //TODO
};

void ElectricMesh::compute_electric_field()
{

    //TODO
};

