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

void ConformalParametrization::compute_dirichlet(SpMat &Dirichlet) {
    Eigen::SparseMatrix<double> C;
    igl::cotmatrix(V, F, C); // Cotangent matrix for the mesh
    int n = V.rows();
    Dirichlet.resize(2 * n, 2 * n);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(2 * C.nonZeros());

    for (int i = 0; i < C.outerSize(); ++i)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(C, i); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
            triplets.emplace_back(it.row() + n, it.col() + n, it.value());
        }
    }
    Dirichlet.setFromTriplets(triplets.begin(), triplets.end());
    Dirichlet.makeCompressed();
}


void ConformalParametrization::compute_area(SpMat &Area)
{
    int n = V.rows();
    Area.resize(2 * n, 2 * n);
    std::vector<Eigen::Triplet<double>> triplets;

    for (int f = 0; f < F.rows(); ++f)
    {
        int i = F(f, 0);
        int j = F(f, 1);
        int k = F(f, 2);

        int idx[3] = {i, j, k};
        for (int e = 0; e < 3; ++e)
        {
            int a = idx[e];
            int b = idx[(e + 1) % 3];

            triplets.emplace_back(a, b + n, 0.5);
            triplets.emplace_back(b, a + n, -0.5);
        }
    }
    Area.setFromTriplets(triplets.begin(), triplets.end());
    Area.makeCompressed();
}

void ConformalParametrization::compute_conformal_energy(SpMat &ConformalEnergy)
{
    ConformalEnergy = Dirichlet - Area;
}



void ConformalParametrization::minimize_energy_spectral(Eigen::MatrixXd &V_uv) {

    int n = V.rows();

    B.resize(2*n, 2*n);
    B.setZero();
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd e = Eigen::VectorXd::Zero(2*n);

    for (int i = 0; i < boundary.size(); ++i) {
        int v = boundary(i);
        e(v) = 1;
        e(v + n) = 1;
        triplets.emplace_back(v, v, 1);
        triplets.emplace_back(v+n, v+n, 1);
    }
    B.setFromTriplets(triplets.begin(), triplets.end());
    int dP = boundary.size();
    Eigen::MatrixXd G = B - (1.0 / dP) * (B * e * e.transpose() * B.transpose());
    //Eigen::MatrixXd G = B;

    Eigen::VectorXd u = Eigen::VectorXd::Random(2 * n);
    u.normalize();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> cholesky_solver;
    cholesky_solver.compute(ConformalEnergy);
    Eigen::VectorXd residual(2 * n);
    int max_iterations = 1000;
    int iterations = 0;
    double tol = 1e-6;
    while (iterations < max_iterations) {
        Eigen::VectorXd w = G * u;
        Eigen::VectorXd v = cholesky_solver.solve(w);
        u = v / v.norm();
        residual = ConformalEnergy * u - (u.dot(ConformalEnergy * u)/u.dot(G * u)) * G * u;
        if (residual.norm() < tol)
            break;
        iterations = iterations + 1;
    }
    V_uv.resize(n, 2);
    V_uv.col(0) = u.segment(0, n);
    V_uv.col(1) = u.segment(n, n);
    Eigen::VectorXd min_uv = V_uv.colwise().minCoeff();
    Eigen::VectorXd max_uv = V_uv.colwise().maxCoeff();
    V_uv = (V_uv.rowwise() - min_uv.transpose()).array().rowwise() / (max_uv - min_uv).transpose().array();
}



void ConformalParametrization::build_parametrizations(){
	compute_dirichlet(Dirichlet);

	compute_area(Area);

	compute_conformal_energy(ConformalEnergy);

	minimize_energy_spectral(V_uv);
}




