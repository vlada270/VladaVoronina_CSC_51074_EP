#include "Eigen/Core"
using Eigen::MatrixXd;
#include <Eigen/Dense>
class CubicInterpolation
{
    MatrixXd M;// coefficients
    MatrixXd a;// solution
    MatrixXd V;// points to interpolate
    MatrixXd y;// left size of linear system

    /*
    Initialize system constraints
    */
    void init_system()
    {
        int n = V.rows();
        for (int i = 0; i < n - 1; ++i)
        {
            double x0 = V(i,0);
            double x1 = V(i+1,0);
            double y0 = V(i,1);
            double y1 = V(i+1,1);

            M(4 * i, 4 * i + 0) = 1;
            M(4 * i, 4 * i + 1) = x0;
            M(4 * i, 4 * i + 2) = x0 * x0;
            M(4 * i, 4 * i + 3) = x0 * x0 * x0;
            y(4 * i) = y0;

            M(4 * i + 1, 4 * i + 0) = 1;
            M(4 * i + 1, 4 * i + 1) = x1;
            M(4 * i + 1, 4 * i + 2) = x1 * x1;
            M(4 * i + 1, 4 * i + 3) = x1 * x1 * x1;
            y(4 * i + 1) = y1;

            if (i > 0) {
                M(4 * i + 2, 4 * (i - 1) + 1) = 1;
                M(4 * i + 2, 4 * (i - 1) + 2) = 2 * x1;
                M(4 * i + 2, 4 * (i - 1) + 3) = 3 * x1 * x1;
                M(4 * i + 2, 4 * i + 1) = -1;
                M(4 * i + 2, 4 * i + 2) = -2 * x1;
                M(4 * i + 2, 4 * i + 3) = -3 * x1 * x1;

                M(4 * i + 3, 4 * (i - 1) + 2) = 2;
                M(4 * i + 3, 4 * (i - 1) + 3) = 6 * x1;
                M(4 * i + 3, 4 * i + 2) = -2;
                M(4 * i + 3, 4 * i + 3) = -6 * x1;
            }
        }
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
        a = M.colPivHouseholderQr().solve(y);

    }


    /*
    Evaluate tangent at step i
    */
    void eval_tangent(int i, MatrixXd &dX, float x)
    {
        // Coefficients for the cubic polynomial in the i-th segment
        double b = a(4 * i + 1);
        double c = a(4 * i + 2);
        double d = a(4 * i + 3);

        // Compute the derivative of the cubic spline (f'(x) = b + 2*c*x + 3*d*x^2)
        dX(i, 0) = 1;
        dX(i, 1) = b + 2 * c * x + 3 * d * x * x;

    }

    /*
    Evaluate function at time t
    */
    double eval_function(double t) {

        int n = V.rows();
        if (t > V(n-1, 0))
            return V(n-1, 1);
        for (int i = 0; i < n - 1; ++i)
        {
            double x0 = V(i, 0);
            double x1 = V(i + 1, 0);

            if (t >= x0 && t <= x1) {
                double result = a(4 * i) + a(4 * i + 1) * t + a(4 * i + 2) * t * t + a(4 * i + 3) * t * t * t;
                return result;
            }
        }
        return 0;
    }
};



