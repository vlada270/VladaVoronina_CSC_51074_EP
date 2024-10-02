
#ifndef UTIL_H
#define UTIL_H
#include "Eigen/Core"

namespace util{
  using namespace Eigen;
  void build_linspace(MatrixXd &linspace, const MatrixXd &V)
  {
    std::cout<<"entering build_linspace"<<std::endl;
    for (size_t i = 0; i < linspace.rows(); i++)
    {
      linspace(i, 0) = V.col(0).minCoeff() + ((V.col(0).maxCoeff() - V.col(0).minCoeff()) / (linspace.rows() - 1)) * i;
    }
    std::cout<<"exiting build_linspace"<<std::endl;
  }



  void add_curve(const MatrixXd &V, std::vector<Vector3d>& points, std::vector<std::array<int, 2>>& edges,std::vector<std::array<double, 3>>& colors, std::array<double, 3> color = {0.0, 0.0, 0.0}){
    // Convert points to a suitable format for Polyscope
    int n_pts = points.size();
    for(int i = 0; i < V.rows(); i++){
      points.push_back(Eigen::Vector3d(V(i, 0), V(i, 1), V(i, 2)));
      if(i<V.rows()-1){
        edges.push_back({n_pts+i, n_pts + (i+1)});
        colors.push_back(color);
      }
    }
  }
  

  
  void translate_points(MatrixXd &V, MatrixXd a){
    for(int i = 0; i < V.rows(); i++){
      V.row(i) += a;
    }
  }
  
}

#endif