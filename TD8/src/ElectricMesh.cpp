#include "ElectricMesh.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/per_vertex_normals.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace Eigen;

ElectricMesh::ElectricMesh(MatrixXd V, MatrixXi F) : Mesh(V,F) {
    set_frame_field();
}

void ElectricMesh::set_frame_field() {
    // Calculate per-vertex normals
    MatrixXd N;
    igl::per_vertex_normals(V, F, N);

    // Initialize the basis vectors for each face
    basisX.resize(F.rows());
    basisY.resize(F.rows());
    bases_tangent_spaces.resize(F.rows());
    bases_tangent_spaces_mat.resize(F.rows());

    for(int i = 0; i < F.rows(); i++) {
        // Get face vertices
        Vector3d v0 = V.row(F(i,0));
        Vector3d v1 = V.row(F(i,1));
        Vector3d v2 = V.row(F(i,2));

        // Compute face normal
        Vector3d normal = (v1 - v0).cross(v2 - v0).normalized();

        // Compute first basis vector (along one edge)
        basisX[i] = (v1 - v0).normalized();

        // Compute second basis vector (perpendicular to both normal and basisX)
        basisY[i] = normal.cross(basisX[i]).normalized();

        // Store bases in matrix form
        MatrixXd basis_mat(2,3);
        basis_mat.row(0) = basisX[i];
        basis_mat.row(1) = basisY[i];
        bases_tangent_spaces_mat[i] = basis_mat;

        // Store as pair of vectors
        bases_tangent_spaces[i] = std::make_pair(basisX[i], basisY[i]);
    }
}

void ElectricMesh::initialize_charge_density(const MatrixXd rho_in) {
    // Store the input charge density
    rho = rho_in;

    // Convert matrix to vector form if needed
    rho_vector.resize(rho.rows());
    for(int i = 0; i < rho.rows(); i++) {
        rho_vector[i] = rho(i,0);
    }
}

void ElectricMesh::solve_for_u() {
    // Compute the cotangent Laplacian matrix
    SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    // Compute the mass matrix
    SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

    // Solve the Poisson equation: -Lu = Mρ
    // Convert to dense format for direct solving
    MatrixXd rho_dense = Map<MatrixXd>(rho_vector.data(), rho_vector.size(), 1);
    MatrixXd b = M * rho_dense;

    // Solve the system
    SimplicialLDLT<SparseMatrix<double>> solver;
    solver.compute(-L);
    u = solver.solve(b);
}

void ElectricMesh::compute_electric_field() {
    // Initialize electric field storage
    electric_field = MatrixXd::Zero(F.rows(), 3);
    electric_field_in_frame.resize(F.rows());

    // For each face
    for(int i = 0; i < F.rows(); i++) {
        // Get vertices of the face
        Vector3d v0 = V.row(F(i,0));
        Vector3d v1 = V.row(F(i,1));
        Vector3d v2 = V.row(F(i,2));

        // Get potential values at vertices
        double u0 = u(F(i,0));
        double u1 = u(F(i,1));
        double u2 = u(F(i,2));

        // Compute gradient in 3D
        Vector3d e1 = v1 - v0;
        Vector3d e2 = v2 - v0;
        double area = e1.cross(e2).norm() / 2.0;

        // Compute electric field (negative gradient of potential)
        Vector3d grad_u = -(u1 - u0) * e1 + (u2 - u0) * e2;
        electric_field.row(i) = grad_u / (2.0 * area);

        // Project electric field onto local frame
        Vector2d E_local;
        E_local(0) = electric_field.row(i).dot(basisX[i]);
        E_local(1) = electric_field.row(i).dot(basisY[i]);
        electric_field_in_frame[i] = E_local;
    }
}

