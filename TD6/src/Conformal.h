#ifndef CONFORMAL_H
#define CONFORMAL_H

#include "Eigen/Core"
#include "Eigen/Sparse"




typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double


class ConformalParametrization {

	private:
		Eigen::VectorXi boundary_flat;

		double residual(const Eigen::VectorXd&);
		// void step_inverse_power_iteration(const Eigen::SimplicialLDLT<SpMat>&, Eigen::VectorXd& );
	public:
        

		Eigen::MatrixXd V_uv,V_uv_spectral,V_uv_harmonic, V,bc;

		Eigen::MatrixXi F;
		SpMat Dirichlet, Area, ConformalEnergy, B, E_b;
		
		Eigen::VectorXi boundary,b;

		ConformalParametrization(const Eigen::MatrixXd &V0, const Eigen::MatrixXi &F0);
		
        double cotangent(double x){
			return 1/(tan(x));
		}

		//compute the matrix for dirichlet energy of the system
		void compute_dirichlet(SpMat &Dirichlet);

		void compute_area(SpMat &Area);

		void compute_conformal_energy(SpMat &ConformalEnergy);

		void minimize_energy_spectral(Eigen::MatrixXd &V_uv);

		void build_parametrizations();

};





#endif // CONFORMAL_H