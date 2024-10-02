#ifndef ELECTRICMESH_H
#define ELECTRICMESH_H
#include "Mesh.h"



class ElectricMesh : public Mesh
{

    private:
        void set_frame_field(); 
        


    public:
        ElectricMesh(Eigen::MatrixXd V, Eigen::MatrixXi F);

        Eigen::MatrixXd rho; // charge density
        std::vector<double> rho_vector;
        Eigen::MatrixXd u;   //

        Eigen::MatrixXd electric_field;
        std::vector<Eigen::Vector2d> electric_field_in_frame;

        std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>> bases_tangent_spaces; // contains for every face a tuple of vectors forming a basis of the tangent space
        std::vector<Eigen::MatrixXd> bases_tangent_spaces_mat; // contains for every face a matrix encoding the basis of the tangent space. The vectors spanning the face are 3D. Then the stored matrix is of size 2 x 3
        std::vector<Eigen::Vector3d> basisX;
        std::vector<Eigen::Vector3d> basisY;
        void initialize_charge_density(const Eigen::MatrixXd rho_in);
        void solve_for_u();
        void compute_electric_field();

        
};

#endif