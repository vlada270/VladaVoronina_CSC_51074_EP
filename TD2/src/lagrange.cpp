#include "Eigen/Dense"
using Eigen::MatrixXd;

class LagrangeInterpolation{
  MatrixXd V;
  MatrixXd W;
  MatrixXd a;

  void Vandermonde(){
      //Complete
      int m = V.rows();
      W = MatrixXd::Ones(m, m);

      for (int i = 0; i < m; ++i) {
          for (int j = 1; j < m; ++j) {
              W(i, j) = pow(V(i, 0), j);  // x^j
          }
      }
      a = W.colPivHouseholderQr().solve(V.col(1));
  }

public:
  LagrangeInterpolation() = default;
  ~LagrangeInterpolation() = default;

  LagrangeInterpolation(const MatrixXd &V1){
    //TODO
      V = V1;
      Vandermonde();
  }

  // evaluate interpolated function at time t
  double eval_function(double t){
      float total = 0;
      //TODO
      int m = a.size();
      for (int i = 0; i < m; i++) {
          total += a(i) * pow(t, i);
      }
      return total;
  }
};
