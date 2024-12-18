#include "gtest/gtest.h"
#include "Eigen/Core"
#include "../src/ElectricMesh.h"
#include <ostream>

#include "igl/readOBJ.h"
#include "igl/readOFF.h"
#include "igl/gaussian_curvature.h"


Eigen::MatrixXd V0;
Eigen::MatrixXi F0;



TEST(CHECK_GRADIENT,ElectricTest){
	#ifdef USING_GNU
	// std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
	int i0 = igl::readOBJ("../data/bunny_new.off", V0, F0);
// Load GNU-specific meshes
#elif defined(USING_MSVC)
	int i0 = igl::readOBJ("../../../../data/bunny_new.off", V0, F0);
	// Load MSVC-specific meshes
#else
	// std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
	int i0 = igl::readOBJ("../data/bunny_new.off", V0, F0);
	// Load default meshes
#endif
	ElectricMesh test_electric_mesh(V0,F0);
	Eigen::MatrixXd rho_in = Eigen::MatrixXd::Zero(V0.rows(), 1);
	rho_in(97, 0) = 1.;
	rho_in(514,0) = -1;
	rho_in(266, 0) = -1.;
	rho_in(110, 0) = 1.;
	test_electric_mesh.initialize_charge_density(rho_in);
	test_electric_mesh.solve_for_u();
	test_electric_mesh.compute_electric_field();
	Eigen::MatrixXd electric_field = test_electric_mesh.electric_field;
	Eigen::MatrixXd test_vector = electric_field.row(200);
	// -0.926459  0.371522 0.0603703
	ASSERT_NEAR(test_vector(0),-0.926459,1e-6);
	ASSERT_NEAR(test_vector(1),0.371522,1e-6);
	ASSERT_NEAR(test_vector(2),0.0603703,1e-6);

}
//Go ahead with the tests

TEST(CHECK_GRADIENT,ElectricTest2){
	#ifdef USING_GNU
		// std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
		int i0 = igl::readOBJ("../data/spot.obj", V0, F0);
	// Load GNU-specific meshes
	#elif defined(USING_MSVC)
		int i0 = igl::readOBJ("../../../../data/spot.obj", V0, F0);
		// Load MSVC-specific meshes
	#else
		// std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
		int i0 = igl::readOBJ("../data/spot.obj", V0, F0);
		// Load default meshes
	#endif

	ElectricMesh test_electric_mesh(V0,F0);
	Eigen::MatrixXd rho_in = Eigen::MatrixXd::Zero(V0.rows(), 1);
	rho_in(97, 0) = 1.;
	rho_in(514,0) = -1;
	rho_in(266, 0) = -1.;
	rho_in(110, 0) = 1.;
	test_electric_mesh.initialize_charge_density(rho_in);
	test_electric_mesh.solve_for_u();
	test_electric_mesh.compute_electric_field();

	double u_check = test_electric_mesh.u(27);
	ASSERT_NEAR(u_check,0.000124508,1e-6);

}
