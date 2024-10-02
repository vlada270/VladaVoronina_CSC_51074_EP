#include "Eigen/Core"

#include <vector>


/**
 * A class for dealin with the computation of Bezier curves and their properties
 **/
class Bezier
{

public:
  /**
   * An iterative implementation of the De Casteljau algorithm for computing a Bezier curve
   *
   * @param V  the vertices of the control polygon
   * @param t   an input parameter in [0, 1]
   *
   * @return    the point B(t), obtaining by evaluating the curve for the value 't'
   **/
  Eigen::MatrixXd de_casteljau(const Eigen::MatrixXd &V, double t);

  Eigen::MatrixXd de_casteljau_intermediate(const Eigen::MatrixXd &V, double t);

  /**
   * Plot the curve, for t=0, dt, 2*dt, ..., 1, with a given 'resolution' <br>
   * where dt=1/(resolution-1)
   *
   * @param resolution  number of points to be evaluated on the curve
   */
  Eigen::MatrixXd plot_curve(const Eigen::MatrixXd &V, int resolution);

  Eigen::MatrixXd plot_curve(const Eigen::MatrixXd &V, int resolution, double t);


  /**
   * Perform the subdivision (once) of the Bezier curve (with parameter t) <br>
   * Return two Bezier curves (with 'n' control points each)
   */
  std::vector<Eigen::MatrixXd> subdivide(const Eigen::MatrixXd &V, double t);
  

  /**
   * Plot the curve using recursive subdivision <br>
   *
   * @param levels  number of levels of subdivisions
   * @return  a polyline representing the curve to be rendered: this is obtained by concantenation of
   * the control points of all subdivided curves
   */
  Eigen::MatrixXd subdivision_plot(const Eigen::MatrixXd &V, int levels);

  /**
   * Compute the tangent of a given curve c(t) for a given parameter t0
   *
   * @param V  the vertices of the control polygon
   * @param t0   an input parameter in [0, 1]
   *
   * @return    the tangent at c(t0)
   **/
  Eigen::MatrixXd compute_tangent(const Eigen::MatrixXd &V, double t0);

  /**
   * Compute the normal vector of a given curve c(t) for a given parameter t0
   *
   * @param V  the vertices of the control polygon
   * @param t0   an input parameter in [0, 1]
   *
   * @return    the normal at c(t0)
   **/
  Eigen::MatrixXd compute_normal(const Eigen::MatrixXd &V, double t0);

  /**
   * Compute a loop of points around a curve c(t) for a given parameter t0
   * The points belongs on a circle lying the hyperplane passing through c(t0) and orthogonal to tangent(t0)
   *
   * @param V  the vertices of the control polygon
   * @param t0   an input parameter in [0, 1]
   *
   * @return    a loop of vertices on the hyperplane passing through c(t0) and orthogonal to tangent(t0)
   **/
  Eigen::MatrixXd compute_loop_of_vertices(const Eigen::MatrixXd &V, double t0, int k, double radius);
  
  //helper function to facilitate the computations
  Eigen::MatrixXd toMatrix(Eigen::Vector3d v){Eigen::MatrixXd M = (Eigen::MatrixXd(1,3)<<v[0],v[1],v[2]).finished();return M;}
  //helper function to facilitate the computations
  Eigen::Vector3d toVector(Eigen::MatrixXd M){Eigen::Vector3d v = Eigen::Vector3d( M(0,0),M(0,1),M(0,2));return v;}
};
