#include "Eigen/Core"

using Eigen::MatrixXd;

class LinearInterpolation{
  MatrixXd V;
  public:
    LinearInterpolation() = default;
    ~LinearInterpolation() = default;

    LinearInterpolation(const MatrixXd &V0){
      V = V0;
    }

    // evaluate function at time t
    double eval_function(double t){
      //complete here

      return 0.;
    }
  };
