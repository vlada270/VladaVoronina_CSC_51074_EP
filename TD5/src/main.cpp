// std includes
#include <iostream>
#include <ostream>
#include <memory>
#include <string>

// libigl includes
#include "igl/readOFF.h"
#include "igl/writeOFF.h"
#include "igl/readOBJ.h"
#include "igl/writeOBJ.h"

// Eigen includes
#include "Eigen/Core"

//Polyscope includes
#include "polyscope/polyscope.h"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"


#include "pca.h"
#include "ICP.h"

//using namespace Eigen; // to use the classes provided by Eigen library
Eigen::MatrixXd V1; // matrix storing vertex coordinates of the input curve
Eigen::MatrixXi F1;
Eigen::MatrixXd V2; // matrix storing vertex coordinates of the input curve
Eigen::MatrixXi F2;

Eigen::MatrixXd normals;

polyscope::PointCloud* psCloud;
polyscope::SurfaceMesh* mesh1;
polyscope::SurfaceMesh* mesh2;

int step =1;
int exercise = 1;


void ex1() {
	#ifdef USING_GNU
		igl::readOFF("../data/star.off", V1, F1);
		igl::readOFF("../data/star_rotated.off", V2, F2);
		
		// Load GNU-specific meshes
	#elif defined(USING_MSVC)
		igl::readOFF("../../../../data/star.off", V1, F1);
		igl::readOFF("../../../../data/star_rotated.off", V2, F2);
		// Load MSVC-specific meshes
	#else
		igl::readOFF("../data/star.off", V1, F1);
		igl::readOFF("../data/star_rotated.off", V2, F2);

		// Load default meshes
	#endif


}

void ex2() {

	#ifdef USING_GNU
		igl::readOFF("../data/bunny.off", V1, F1);
		igl::readOFF("../data/bunny_rotated.off", V2, F2);
		// Load GNU-specific meshes
	#elif defined(USING_MSVC)
		igl::readOFF("../../../../data/bunny.off", V1, F1);
		igl::readOFF("../../../../data/bunny_rotated.off", V2, F2);
		// Load MSVC-specific meshes
	#else
		igl::readOFF("../data/bunny.off", V1, F1);
		igl::readOFF("../data/bunny_rotated.off", V2, F2);
		// Load default meshes
	#endif
}

void ex3() {
	#ifdef USING_GNU
		igl::readOFF("../data/bunny.off", V1, F1);
		igl::readOFF("../data/bunny_rotated.off", V2, F2);
		// Load GNU-specific meshes
	#elif defined(USING_MSVC)
		igl::readOFF("../../../../data/bunny.off", V1, F1);
		igl::readOFF("../../../../data/bunny_rotated.off", V2, F2);
		// Load MSVC-specific meshes
	#else
		igl::readOFF("../data/bunny.off", V1, F1);
		igl::readOFF("../data/bunny_rotated.off", V2, F2);
		// Load default meshes
	#endif

}

void ex4() {
	#ifdef USING_GNU
		igl::readOFF("../data/sphere.off", V1, F1);
		// Load GNU-specific meshes
	#elif defined(USING_MSVC)
		igl::readOFF("../../../../data/sphere.off", V1, F1);
		// Load MSVC-specific meshes
	#else
		igl::readOFF("../data/sphere.off", V1, F1);
	#endif

}

void callback() {
	static int numPoints = 2000;
	static float param = 3.14;

	
	if(ImGui::SliderInt("Exercise", &step, 1, 4)) {}

	if (ImGui::Button("Load Mesh")) {
		if (step == 1) {
			ex1();
			mesh1 = polyscope::registerSurfaceMesh("ICP Mesh 1", V1, F1);
			mesh2 = polyscope::registerSurfaceMesh("ICP Mesh 2", V2, F2);
		}
		else if (step == 2) {
			ex2();
			mesh1 = polyscope::registerSurfaceMesh("ICP Mesh 1", V1, F1);
			mesh2 = polyscope::registerSurfaceMesh("ICP Mesh 2", V2, F2);
		}
		else if (step == 3) {
			ex3();
			mesh1 = polyscope::registerSurfaceMesh("ICP Mesh 1", V1, F1);
			mesh2 = polyscope::registerSurfaceMesh("ICP Mesh 2", V2, F2);
		}
		else if (step == 4) {
			ex4();
			psCloud = polyscope::registerPointCloud("Point Cloud", V1);
		}
	}	
		

	ImGui::PushItemWidth(100);
	if (ImGui::Button("ICP Step Ex1")) {
		ICP::transform(V1,V2);
		polyscope::registerSurfaceMesh("ICP Mesh 1", V1, F1);
	}

	if (ImGui::Button("ICP Ex2")) {
    
		Eigen::MatrixXd nn_V2(V1.rows(), V1.cols());
		ICP::nearest_neighbour(V1, V2, nn_V2, ICP::KnnStrategy::OCTREE);
		ICP::transform(V1, nn_V2);
		polyscope::registerSurfaceMesh("ICP Mesh 1", V1, F1);
      
	}
	
	if (ImGui::Button("ICP Point to Plane Ex3")) {
		Eigen::MatrixXd nn_V2(V1.rows(), V1.cols());    
		ICP::nearest_neighbour_point_to_plane(V1, V2, nn_V2);
		//  complete here by displaying the pair wise distances
		ICP::transform(V1,nn_V2);
		polyscope::registerSurfaceMesh("ICP Mesh 1", V1, F1);
	}
	
	if (ImGui::Button("PCA Ex4")) {
		
		Eigen::MatrixXi I;
		PCA::k_nearest_neighbour(V1,I,12);
		PCA::compute_normals(V1,I, 12, normals);
		psCloud->addVectorQuantity("normals", normals);
	}

	ImGui::PopItemWidth();
}






int main(int argc, char *argv[])
{
	std::cout<<" -- Iterative Closest Points and Principle Component Analysis -- "<<std::endl;

  /*if (argc != 2) {
    std::cerr << "./td5_executable <mode>, mode = 1/2/3/4" << std::endl;
    return 1;
  }
  if (!strcmp(argv[1],"1")) {
    ex1();
  }
  if (!strcmp(argv[1],"2")) {
    ex2();
  }
  if (!strcmp(argv[1],"3")) {
    ex3();
  }
  if (!strcmp(argv[1],"4")) {
    ex4();
  }*/
	// --- Visualize the mesh using Polyscope ---
	//--------------------------------------------

	// Options
	polyscope::options::autocenterStructures = true;
	polyscope::view::windowWidth = 1024;
	polyscope::view::windowHeight = 1024;

	// Initialize polyscope
	polyscope::init();

  
	 
	

	polyscope::state::userCallback = callback;


	// Show the gui
	polyscope::show();

	return 0;

}
