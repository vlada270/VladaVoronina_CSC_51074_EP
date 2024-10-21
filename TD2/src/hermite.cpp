#include "Eigen/Core"
#include "Eigen/LU"

using Eigen::MatrixXd;

class HermiteInterpolation {

    MatrixXd solutionx; // solution for x
    MatrixXd solutiony; // solution for y
    MatrixXd steps;     // points to interpolate
    MatrixXd solvey;    // right side of linear system for y
    MatrixXd solvex;    // right side of linear system for x
    MatrixXd slopes;
    MatrixXd W;

    /*
    Solve for the x(t) cubic polynomial coefficients
    */
    void solve_x() {
    }

    /*
    Solve for the y(t) cubic polynomial coefficients
    */
    void solve_y() {
    }

public:
    HermiteInterpolation() = default;
    ~HermiteInterpolation() = default;


    HermiteInterpolation(const MatrixXd &V1) {
        // complete here

    }

    /*
    Evaluate function at parameter t (generate interpolated points)
    */
    void eval_function(MatrixXd &linspace) {
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

