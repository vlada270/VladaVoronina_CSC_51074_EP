

#define EIGEN_NO_STATIC_ASSERT
#include "ICP.h"
#include "pca.h"
#include <thread>
#include <chrono>
#include "Eigen/SVD"
#include "Eigen/Eigenvalues"
using namespace Eigen;
#include "igl/octree.h"
#include "igl/knn.h"
#include <igl/per_vertex_normals.h>


void ICP::nearest_neighbour(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2, KnnStrategy strategy) {
    if (strategy == KnnStrategy::OCTREE) {
        //TODO
        std::vector<std::vector<int>> point_indices;
        Eigen::MatrixXi CH;
        Eigen::MatrixXd CN;
        Eigen::VectorXd W;

        igl::octree(V2, point_indices, CH, CN, W);
        Eigen::MatrixXi closest_indices;
        igl::knn(V1, V2, 1, point_indices, CH, CN, W, closest_indices);
        nn_V2.resize(V1.rows(), V1.cols());
        for (int i = 0; i < V1.rows(); ++i) {
            nn_V2.row(i) = V2.row(closest_indices(i, 0));
        }
    } else {
        //TODO
        nn_V2.resize(V1.rows(), V1.cols());
        for (int i = 0; i < V1.rows(); ++i) {
            double min_dist = std::numeric_limits<double>::max();
            int min_index = -1;
            for (int j = 0; j < V2.rows(); ++j) {
                double dist = (V1.row(i) - V2.row(j)).squaredNorm();
                if (dist < min_dist) {
                    min_dist = dist;
                    min_index = j;
                }
            }
            nn_V2.row(i) = V2.row(min_index);
        }
    }
}

void ICP::nearest_neighbour_point_to_plane(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2) {
    // return the nearest neighbour to V1 in V2 as nn_V2 using the point to plane algorithm

    //TODO
    nn_V2.resize(V1.rows(), V1.cols());

    std::vector<std::vector<int>> point_indices;
    Eigen::MatrixXi CH;
    Eigen::MatrixXd CN;
    Eigen::VectorXd W;
    igl::octree(V2, point_indices, CH, CN, W);
    MatrixXd V2_normals(V2.rows(), 3);
    MatrixXi closest_indices;
    igl::knn(V2, V2, 10, point_indices, CH, CN, W, closest_indices);

    for (int i = 0; i < V2.rows(); i++) {
        MatrixXd neighbors(10, 3);
        for (int j = 0; j < 10; j++) {
            neighbors.row(j) = V2.row(closest_indices(i, j));
        }

        Vector3d centroid = neighbors.colwise().mean();
        MatrixXd centered = neighbors.rowwise() - centroid.transpose();
        Matrix3d covariance = centered.transpose() * centered;
        SelfAdjointEigenSolver<Matrix3d> eig(covariance);
        V2_normals.row(i) = eig.eigenvectors().col(0).normalized();
    }

    igl::knn(V1, V2, 1, point_indices, CH, CN, W, closest_indices);

    for (int i = 0; i < V1.rows(); ++i) {
        Vector3d p = V1.row(i);
        int closest_idx = closest_indices(i, 0);
        Vector3d q = V2.row(closest_idx);
        Vector3d n = V2_normals.row(closest_idx);
        nn_V2.row(i) = V2.row(closest_idx);
    }
}

void ICP::transform(MatrixXd &V1, const MatrixXd &V2) {
    //align V1 to V2 when V1 and V2 points are in correspondance
    //TODO
    Eigen::Vector3d centroid_V1 = V1.colwise().mean();
    Eigen::Vector3d centroid_V2 = V2.colwise().mean();
    Eigen::MatrixXd centered_V1 = V1.rowwise() - centroid_V1.transpose();
    Eigen::MatrixXd centered_V2 = V2.rowwise() - centroid_V2.transpose();
    Eigen::MatrixXd H = centered_V1.transpose() * centered_V2;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::Matrix3d R = V * U.transpose();

    if (R.determinant() < 0) {
        V.col(2) *= -1;
        R = V * U.transpose();
    }

    Eigen::Vector3d t = centroid_V2 - R * centroid_V1;
    V1 = (R * V1.transpose()).colwise() + t;
    V1.transposeInPlace();
}

