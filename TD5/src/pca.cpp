#include "pca.h"
#include "Eigen/SVD"
#include "Eigen/Eigenvalues"
#include "igl/octree.h"
#include "igl/knn.h"

using namespace Eigen;

void PCA::k_nearest_neighbour(const MatrixXd &V1,Eigen::MatrixXi &I, int k){
  // Build octree

}

void PCA::compute_normals(const MatrixXd &V1,const Eigen::MatrixXi &I, int k, MatrixXd &normals){
  //TODO


}





