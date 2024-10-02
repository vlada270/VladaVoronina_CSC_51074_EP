#include "igl/boundary_loop.h"
#include "igl/boundary_facets.h"
#include "igl/cotmatrix.h"
#include <Eigen/IterativeLinearSolvers>




#include "Conformal.h"


using Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd, Eigen::VectorXi, Eigen::Vector3d, Eigen::SparseMatrix, Eigen::Triplet, Eigen::Map;


ConformalParametrization::ConformalParametrization(const MatrixXd &V0, const MatrixXi &F0){
	V =  V0;
	F = F0;
	igl::boundary_loop(F,boundary);

	//fix points on the boundary
	b.resize(2,1);
	b(0) = boundary(0);
	b(1) = boundary(boundary.size()/2);
}
		
//compute the matrix for dirichlet energy of the system
void ConformalParametrization::compute_dirichlet(SpMat &Dirichlet){
	//TODO
	
}



void ConformalParametrization::compute_area(SpMat &Area){
	//TODO
	
}

void ConformalParametrization::compute_conformal_energy(SpMat &ConformalEnergy){
	//TODO
	
	
}



void ConformalParametrization::minimize_energy_spectral(Eigen::MatrixXd &V_uv){
	
	//TODO

}

void ConformalParametrization::build_parametrizations(){
	compute_dirichlet(Dirichlet);

	compute_area(Area);

	compute_conformal_energy(ConformalEnergy);

	minimize_energy_spectral(V_uv);
}




