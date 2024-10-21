#include "Eigen/Core"
#include <iostream>

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
      int n = V.cols();
      int m = V.rows();

      if (t > V(m-1, 0))
          return V(m-1, 1);
      for (int i = 0; i < m - 1; i++) {
          if (t >= V(i, 0) && t <= V(i+1, 0)) {
              double x1 = V(i, 0);
              double x2 = V(i + 1, 0);
              double y1 = V(i, 1);
              double y2 = V(i+1, 1);
              double eval = y1 + (t - x1) * ((y2 - y1) / (x2 - x1));
              return eval;
          }
      }
      return 0.;
    }
  };
