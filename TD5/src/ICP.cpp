

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


void ICP::nearest_neighbour(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2,  ICP::KnnStrategy strategy = ICP::KnnStrategy::OCTREE){
  if (strategy == KnnStrategy::OCTREE) {
    //TODO
  } else {
    //TODO
  }
}

void ICP::nearest_neighbour_point_to_plane(const MatrixXd &V1, const MatrixXd &V2, MatrixXd &nn_V2) {
  // return the nearest neighbour to V1 in V2 as nn_V2 using the point to plane algorithm
	
	//TODO

}

void ICP::transform(MatrixXd &V1,const MatrixXd &V2){
  //align V1 to V2 when V1 and V2 points are in correspondance
  
  //TODO

}
