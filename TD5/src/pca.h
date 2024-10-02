#ifndef PCA_H
#define PCA_H
#include "Eigen/Core"

namespace PCA{
    void k_nearest_neighbour(const Eigen::MatrixXd &V1,Eigen::MatrixXi &I, int k);
    void compute_normals(const Eigen::MatrixXd &V1,const Eigen::MatrixXi &I, int k, Eigen::MatrixXd &normals);
}
#endif