#include "Bezier.h"
#include "Eigen/Geometry"

using Eigen::MatrixXd, Eigen::RowVector3d,Eigen::Vector3d;

using namespace Eigen;

MatrixXd Bezier::de_casteljau(const MatrixXd &V, double t)
{
    MatrixXd p = V;
    for (int i = 1; i < V.rows(); ++i) {
        for (int j = 0; j < V.rows() - i; ++j) {
            p.row(j) = (1 - t) * p.row(j) + t * p.row(j + 1);
        }
    }
    return p.row(0);
}



MatrixXd Bezier::de_casteljau_intermediate(const MatrixXd &V, double t){

    MatrixXd intermediate(V.rows() - 1, V.cols());

    for (int i = 0; i < V.rows() - 1; ++i) {
        intermediate.row(i) = (1 - t) * V.row(i) + t * V.row(i + 1);
    }

    return intermediate;


}

/**
 * Plot the curve, for t=0, dt, 2*dt, ..., 1, with a given 'resolution' <br>
 * where dt=1/(resolution-1)
 *
 * @param resolution  number of points to be evaluated on the curve
 */

MatrixXd Bezier::plot_curve(const MatrixXd &V, int resolution)
{
    MatrixXd m(resolution, 3);
    double dt = 1.0 / (resolution - 1);

    for (int i = 0; i < resolution; ++i) {
        double t = i * dt;
        m.row(i) = de_casteljau(V, t);
    }

    return m;
}

MatrixXd Bezier::plot_curve(const MatrixXd &V, int resolution, double t)
{
    MatrixXd m((V.rows() * (V.rows() - 1)) / 2, V.cols());

    MatrixXd p = V;
    int row = 0;
    for (int i = 1; i < V.rows(); ++i)
    {
        p = de_casteljau_intermediate(p, t);
        for (int j = 0; j < p.rows(); ++j)
            m.row(row++) = p.row(j);
    }
    return m;
}


/**
 * Perform the subdivision (once) of the Bezier curve (with parameter t) <br>
 * Return two Bezier curves (with 'n' control points each)
 */
std::vector<MatrixXd> Bezier::subdivide(const MatrixXd &V, double t)
{
    std::vector<MatrixXd> subDiv(2);
    MatrixXd p = V;
    MatrixXd subdivideLeft(V.rows(), V.cols());
    MatrixXd subdivideRight(V.rows(), V.cols());

    for (int i = 0; i < V.rows(); ++i) {
        subdivideLeft.row(i) = p.row(0);
        subdivideRight.row(V.rows() - 1 - i) = p.row(p.rows() - 1 - i);
        for (int j = 0; j < V.rows() - 1 - i; ++j)
            p.row(j) = (1 - t) * p.row(j) + t * p.row(j + 1);
    }
    subDiv[0] = subdivideLeft;
    subDiv[1] = subdivideRight;
    return subDiv;
}

/**
 * Plot the curve using recursive subdivision <br>
 *
 * @param levels  number of levels of subdivisions
 * @return  a polyline representing the curve to be rendered: this is obtained by concantenation of
 * the control points of all subdivided curves
 */
MatrixXd Bezier::subdivision_plot(const MatrixXd &V, int levels)
{
    if (levels == 0)
        return V;
    std::vector<MatrixXd> subdividedCurves = subdivide(V, 0.5);
    MatrixXd left = subdivision_plot(subdividedCurves[0], levels - 1);
    MatrixXd right = subdivision_plot(subdividedCurves[1], levels - 1);

    MatrixXd result(left.rows() + right.rows(), V.cols());
    result.topRows(left.rows()) = left;
    result.bottomRows(right.rows()) = right;

    return result;
}

/**
 * Compute the tangent of a given curve c(t) for a given parameter t0
 *
 * @param V  the vertices of the control polygon
 * @param t0   an input parameter in [0, 1]
 *
 * @return    the tangent at c(t0)
 **/
MatrixXd Bezier::compute_tangent(const MatrixXd &V, double t0)
{
    MatrixXd derivative(V.rows() - 1, V.cols());
    for (int i = 0; i < V.rows() - 1; ++i)
        derivative.row(i) = (V.rows() - 1) * (V.row(i + 1) - V.row(i));
    return de_casteljau(derivative, t0);
}

/**
 * Compute the normal vector of a given curve c(t) for a given parameter t0
 *
 * @param V  the vertices of the control polygon
 * @param t0   an input parameter in [0, 1]
 *
 * @return    the normal at c(t0)
 **/
MatrixXd Bezier::compute_normal(const MatrixXd &V, double t0)
{
    MatrixXd second(V.rows() - 2, V.cols());
    for (int i = 0; i < V.rows() - 2; ++i)
        second.row(i) = (V.rows() - 2) * (V.rows()-1) * (V.row(i + 2) - 2 * V.row(i + 1) + V.row(i));
    return de_casteljau(second, t0);
}

/**
 * Compute a loop of points around a curve c(t) for a given parameter t0
 * The points belongs on a circle lying the hyperplane passing through c(t0) and orthogonal to tangent(t0)
 *
 * @param V  the vertices of the control polygon
 * @param t0   an input parameter in [0, 1]
 *
 * @return    a loop of vertices on the hyperplane passing through c(t0) and orthogonal to tangent(t0)
 **/
MatrixXd Bezier::compute_loop_of_vertices(const MatrixXd &V, double t0, int k, double radius)
{
    MatrixXd p(k + 1, 3);
    Vector3d circle = de_casteljau(V, t0).transpose();
    Vector3d tangent = compute_tangent(V, t0).transpose().normalized();
    Vector3d orthogonal1, orthogonal2;
    Vector3d x(0, 0, 1);
    double angleIncrement = (2.0 * M_PI) / k;

    if (std::abs(tangent.dot(x)) > 0.95)
        x = Vector3d(0, 1, 0);
    orthogonal1 = tangent.cross(x).normalized();
    orthogonal2 = tangent.cross(orthogonal1).normalized();

    for (int i = 0; i < k; ++i)
        p.row(i) = circle + radius * (std::cos(angleIncrement * i) * orthogonal1 + std::sin(angleIncrement * i) * orthogonal2);
    p.row(k) = p.row(0);
    return p;
}
