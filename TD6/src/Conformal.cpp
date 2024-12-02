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
    int num_vertices = V.rows();
    Dirichlet.resize(2 * num_vertices, 2 * num_vertices);
    std::vector<Eigen::Triplet<double>> triplets;

    for (int f = 0; f < F.rows(); ++f) {
        int v0 = F(f, 0);
        int v1 = F(f, 1);
        int v2 = F(f, 2);

        Eigen::Vector3d p0 = V.row(v0);
        Eigen::Vector3d p1 = V.row(v1);
        Eigen::Vector3d p2 = V.row(v2);

        Eigen::Vector3d e0 = p1 - p0;
        Eigen::Vector3d e1 = p2 - p1;
        Eigen::Vector3d e2 = p0 - p2;

        double angle0 = std::acos(e1.dot(-e2) / (e1.norm() * e2.norm()));
        double angle1 = std::acos(e2.dot(-e0) / (e2.norm() * e0.norm()));
        double angle2 = std::acos(e0.dot(-e1) / (e0.norm() * e1.norm()));

        double cot0 = cotangent(angle0);
        double cot1 = cotangent(angle1);
        double cot2 = cotangent(angle2);


        triplets.push_back(Eigen::Triplet<double>(v0, v0, cot1 + cot2));
        triplets.push_back(Eigen::Triplet<double>(v1, v1, cot2 + cot0));
        triplets.push_back(Eigen::Triplet<double>(v2, v2, cot0 + cot1));

        triplets.push_back(Eigen::Triplet<double>(v0, v1, -cot2));
        triplets.push_back(Eigen::Triplet<double>(v1, v0, -cot2));

        triplets.push_back(Eigen::Triplet<double>(v1, v2, -cot0));
        triplets.push_back(Eigen::Triplet<double>(v2, v1, -cot0));

        triplets.push_back(Eigen::Triplet<double>(v2, v0, -cot1));
        triplets.push_back(Eigen::Triplet<double>(v0, v2, -cot1));

        triplets.push_back(Eigen::Triplet<double>(v0 + num_vertices, v0 + num_vertices, cot1 + cot2));
        triplets.push_back(Eigen::Triplet<double>(v1 + num_vertices, v1 + num_vertices, cot2 + cot0));
        triplets.push_back(Eigen::Triplet<double>(v2 + num_vertices, v2 + num_vertices, cot0 + cot1));

        triplets.push_back(Eigen::Triplet<double>(v0 + num_vertices, v1 + num_vertices, -cot2));
        triplets.push_back(Eigen::Triplet<double>(v1 + num_vertices, v0 + num_vertices, -cot2));

        triplets.push_back(Eigen::Triplet<double>(v1 + num_vertices, v2 + num_vertices, -cot0));
        triplets.push_back(Eigen::Triplet<double>(v2 + num_vertices, v1 + num_vertices, -cot0));

        triplets.push_back(Eigen::Triplet<double>(v2 + num_vertices, v0 + num_vertices, -cot1));
        triplets.push_back(Eigen::Triplet<double>(v0 + num_vertices, v2 + num_vertices, -cot1));
    }

    Dirichlet.setFromTriplets(triplets.begin(), triplets.end());
}


void ConformalParametrization::compute_area(SpMat &Area){
	//TODO
    int num_vertices = V.rows();
    Area.resize(2 * num_vertices, 2 * num_vertices);
    std::vector<Eigen::Triplet<double>> triplets;
    std::vector<double> vertex_areas(num_vertices, 0.0);

    for (int i = 0; i < F.rows(); ++i) {
        Eigen::Vector3d a = V.row(F(i, 1)) - V.row(F(i, 0));
        Eigen::Vector3d b = V.row(F(i, 2)) - V.row(F(i, 0));
        double triangle_area = 0.5 * a.cross(b).norm();

        for (int j = 0; j < 3; ++j) {
            vertex_areas[F(i, j)] += triangle_area / 3.0;
        }
    }

    for (int v = 0; v < num_vertices; ++v) {
        triplets.push_back(Eigen::Triplet<double>(v, v, vertex_areas[v]));
        triplets.push_back(Eigen::Triplet<double>(v + num_vertices, v + num_vertices, vertex_areas[v]));
    }

    Area.setFromTriplets(triplets.begin(), triplets.end());

}

void ConformalParametrization::compute_conformal_energy(SpMat &ConformalEnergy) {
    int num_vertices = V.rows();

    ConformalEnergy.resize(2 * num_vertices, 2 * num_vertices);

    std::vector<Eigen::Triplet<double>> triplets;

    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            int v1 = F(i, j);
            int v2 = F(i, (j + 1) % 3);
            int v3 = F(i, (j + 2) % 3);

            Eigen::Vector3d edge1 = V.row(v2) - V.row(v1);
            Eigen::Vector3d edge2 = V.row(v3) - V.row(v1);

            double l1 = edge1.norm();
            double l2 = edge2.norm();

            double energy_contribution = 4.0 / (l1 * l1 + l2 * l2);

            triplets.push_back(Eigen::Triplet<double>(v1, v1, energy_contribution));

            triplets.push_back(Eigen::Triplet<double>(v1 + num_vertices, v1 + num_vertices, energy_contribution));
        }
    }

    ConformalEnergy.setFromTriplets(triplets.begin(), triplets.end());
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




