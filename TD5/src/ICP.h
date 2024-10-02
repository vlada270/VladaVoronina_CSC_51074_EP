#ifndef ICP_H
#define ICP_H
#include <Eigen/Core>

namespace ICP{
  enum KnnStrategy {
    OCTREE,
    BRUTEFORCE
  };

  void nearest_neighbour(const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2, Eigen::MatrixXd &nn_V2, KnnStrategy strategy);
  void nearest_neighbour_point_to_plane(const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2, Eigen::MatrixXd &nn_V2);
  void transform(Eigen::MatrixXd &V1,const Eigen::MatrixXd &V2);
}
#endif