#include "Eigen/Core"

using Eigen::MatrixXd;

class HermiteInterpolation{

  MatrixXd solutionx;// solution for x
  MatrixXd solutiony;// solution for y
  MatrixXd steps;// points to interpolate
  MatrixXd solvey;// left size of linear system
  MatrixXd solvex;// left size of linear system
  MatrixXd slopes;
  MatrixXd W;

  void solve_x(){
    // complete here
    
  }

  void solve_y(){
    // complete here
    
  }

public:

  HermiteInterpolation() = default;
  ~HermiteInterpolation() = default;

  HermiteInterpolation(const MatrixXd &V1){
    // complete here
    
  }

  // complete linspace with corrdinates
  void eval_function(MatrixXd &linspace){
      // complete here
  }

  /*
  Evaluate tangent at step i
  */
  void eval_tangent(int i, MatrixXd& dX)
  {
    // complete here
  }
};
