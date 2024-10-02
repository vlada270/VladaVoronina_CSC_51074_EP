#ifndef LAPLACIANMESH_H
#define LAPLACIANMESH_H
#include "Mesh.h"
#include "Eigen/Sparse"
#include "Eigen/Core"

typedef Eigen::SparseMatrix<double> SpMat;

class LaplacianMesh : public Mesh
{
private:
    /* data */
    void compute_dirichlet();
    void compute_area_matrix();
    void compute_laplacian();

    double cotangent(double x){
		return 1/(tan(x));
	}

public:
    LaplacianMesh(/* args */);
    LaplacianMesh(Eigen::MatrixXd V,Eigen::MatrixXi F);
    ~LaplacianMesh();
    SpMat L; //cotangent
    SpMat A; //normalization for the laplacian
    SpMat Ainv; //inverse of A for performance
    SpMat Delta; //laplacian

    void setHeat(Eigen::VectorXd & u);
    void heat_step_explicit(Eigen::VectorXd & u, double time_step );
    void heat_step_implicit(Eigen::VectorXd & u, double time_step);    
};




#endif // LAPLACIANMESH_H