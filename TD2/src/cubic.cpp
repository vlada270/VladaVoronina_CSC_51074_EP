#include "Eigen/Core"
using Eigen::MatrixXd;

class CubicInterpolation
{
  MatrixXd M;// coefficients
  MatrixXd a;// solution
  MatrixXd V;// points to interpolate
  MatrixXd y;// left size of linear system

  /*
  Initialize system constraints
  */
  void init_system(){
    // complete here
     
  }

public:

  CubicInterpolation() = default;
  ~CubicInterpolation() = default;

  CubicInterpolation(const MatrixXd &V1){
    M = MatrixXd::Zero(4*(V1.rows() - 1),4*(V1.rows() - 1));
    y = MatrixXd::Zero(4*(V1.rows() - 1),1);
    this->V = V1;
    init_system();
    // complete here
    
  }


  /*
  Evaluate tangent at step i
  */
  void eval_tangent(int i, MatrixXd &dX, float x)
  {
    // complete here
      
  }

  /*
  Evaluate function at time t
  */
  double eval_function(double t){
    // complete here
    return 0;
  }
};
