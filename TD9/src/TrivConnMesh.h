#ifndef TRIVCONNMESH_H
#define TRIVCONNMESH_H

#include "Mesh.h"
#include <Eigen/Sparse>


typedef Eigen::SparseMatrix<double> SparseMatrixXd;

class TrivConnMesh : public Mesh {

    private:
        std::vector<bool> in_tree_information;
        std::vector<VertexTree> spanningTreeDual;

    protected:

        Eigen::VectorXd dihedral_angle; // contains for every primal edge the hinge angle.
        Eigen::VectorXd holonomy; //contains per vertex the holonomy around the dual cycles 

        // SparseMatrixXd A, A_check,d0,G;


        Eigen::VectorXd singularities;
        Eigen::VectorXd target_holonomies;

        std::vector<double> index_to_basis;


        
        

        double adjustToRange(double input);

    public:

        TrivConnMesh(Eigen::MatrixXd, Eigen::MatrixXi );
        TrivConnMesh();

        void initialize_edges();

        void set_frame_field();
        
        void compute_hinge_connection();
        
        void compute_angles_defects();
        
        void compute_spanning_tree();
        
        // Build the linear system for the trivial connections
        void build_RHS_for_problem();

        void build_LHS();

        //run this method to test if the 
        void check_RHS_for_consistency(const Eigen::VectorXd&);
        
        // solve for the adjustment angles.
        void solve_for_angles();

        //transport one chosen vector with the given alignment vector
        void compute_transported_vector_field();


        SparseMatrixXd A;

        std::vector<std::pair<double, double>> parallel_transported_field;
        std::vector<Eigen::Vector2d> parallel_transported_field_vec;

        std::vector<std::array<size_t, 2>> edges_spanning_tree;
        
        

        Eigen::VectorXd adjustment_angles; // contains for every primal edge the adjustment angle to turn the hinge connection into the trivial connection
        Eigen::VectorXd adjustment_angles_frame;

        std::vector<EdgePtr> edges;

        std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>> bases_tangent_spaces; // contains for every face a tuple of vectors forming a basis of the tangent space
        std::vector<Eigen::MatrixXd> bases_tangent_spaces_mat; // contains for every face a matrix encoding the basis of the tangent space. The vectors spanning the face are 3D. Then the stored matrix is of size 2 x 3

        std::vector<Eigen::Vector3d> basisX;
        std::vector<Eigen::Vector3d> basisY;
        std::vector<std::array<size_t, 2>> edges_dual_spanning_tree;
};


#endif