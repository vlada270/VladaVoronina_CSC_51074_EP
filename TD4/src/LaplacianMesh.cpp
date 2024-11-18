#include "LaplacianMesh.h"
#include "Eigen/IterativeLinearSolvers"

using Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd, Eigen::Vector3d;


LaplacianMesh::LaplacianMesh(/* args */){}

LaplacianMesh::~LaplacianMesh(){}

void LaplacianMesh::compute_dirichlet(){
    //TODO
    std::vector<Eigen::Triplet<double>> triplets;

    for (const auto& face : primal_faces) {
        auto v0 = face->getVertices()[0];
        auto v1 = face->getVertices()[1];
        auto v2 = face->getVertices()[2];

        double cot_alpha = cotangent((v1->position - v0->position).normalized().dot((v2->position - v0->position).normalized()));
        double cot_beta = cotangent((v2->position - v1->position).normalized().dot((v0->position - v1->position).normalized()));
        double cot_gamma = cotangent((v0->position - v2->position).normalized().dot((v1->position - v2->position).normalized()));
        triplets.emplace_back(v0->index, v1->index, 0.5 * cot_alpha);
        triplets.emplace_back(v1->index, v0->index, 0.5 * cot_alpha);

        triplets.emplace_back(v1->index, v2->index, 0.5 * cot_beta);
        triplets.emplace_back(v2->index, v1->index, 0.5 * cot_beta);

        triplets.emplace_back(v2->index, v0->index, 0.5 * cot_gamma);
        triplets.emplace_back(v0->index, v2->index, 0.5 * cot_gamma);
        triplets.emplace_back(v0->index, v0->index, -0.5 * (cot_alpha + cot_gamma));
        triplets.emplace_back(v1->index, v1->index, -0.5 * (cot_alpha + cot_beta));
        triplets.emplace_back(v2->index, v2->index, -0.5 * (cot_beta + cot_gamma));
    }

    L.resize(V.rows(), V.rows());
    L.setFromTriplets(triplets.begin(), triplets.end());
}

void LaplacianMesh::compute_area_matrix() {
    std::vector<Eigen::Triplet<double>> triplets;

    for (const auto& vertex : primal_vertices) {
        double area_sum = 0.0;

        for (const auto& face : vertex->getOneRingFaces()) {
            auto v0 = face->getVertices()[0]->position;
            auto v1 = face->getVertices()[1]->position;
            auto v2 = face->getVertices()[2]->position;

            // Compute triangle area
            double area = ((v1 - v0).cross(v2 - v0)).norm() / 2.0;
            area_sum += area / 3.0; // 1/3rd contribution to each vertex
        }

        triplets.emplace_back(vertex->index, vertex->index, area_sum);
    }

    A.resize(V.rows(), V.rows());
    A.setFromTriplets(triplets.begin(), triplets.end());
    Ainv = A;
    for (int k = 0; k < A.rows(); ++k) {
        if (A.coeff(k, k) > 1e-10) {
            Ainv.coeffRef(k, k) = 1.0 / A.coeff(k, k);
        }
    }
}


void LaplacianMesh::compute_laplacian() {
    Delta = Ainv * L;
}


LaplacianMesh::LaplacianMesh(Eigen::MatrixXd V,Eigen::MatrixXi F):Mesh(V,F){
    this->V = V;
    this->F = F;

    compute_dirichlet();

    compute_area_matrix();

    compute_laplacian();

}


void LaplacianMesh::setHeat(Eigen::VectorXd & u) {
    u = Eigen::VectorXd::Zero(V.rows());

    for (int i = 0; i < V.rows(); i++) {
        u(i) = V(i, 0);
    }
    std::cout << "Initialized:\n" << u.transpose() << std::endl;
}

void LaplacianMesh::heat_step_explicit(Eigen::VectorXd& u, double time_step) {
    u = u + time_step * (Delta * u);
}

void LaplacianMesh::heat_step_implicit(Eigen::VectorXd& u, double time_step) {
    Eigen::SparseMatrix<double> I(V.rows(), V.rows());
    I.setIdentity();

    Eigen::SparseMatrix<double> system = I - time_step * Delta;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
    solver.compute(system);
    u = solver.solve(u);

}
