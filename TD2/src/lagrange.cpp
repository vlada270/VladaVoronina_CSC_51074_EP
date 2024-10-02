#include "Eigen/Dense"
using Eigen::MatrixXd;

class LagrangeInterpolation{
  MatrixXd V;
  MatrixXd W;
  MatrixXd a;

  void Vandermonde(){
      //Complete
  }

public:
  LagrangeInterpolation() = default;
  ~LagrangeInterpolation() = default;

  LagrangeInterpolation(const MatrixXd &V1){
    //TODO
  }

  // evaluate interpolated function at time t
  double eval_function(double t){
      float total = 0;
      //TODO
      return total;
  }
};
