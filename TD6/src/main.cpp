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
#include "igl/lscm.h"
#include "igl/boundary_loop.h"


//Polyscope includes
#include "polyscope/polyscope.h"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"



#include "Conformal.h"



Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd V_uv_libigl;
Eigen::MatrixXd V_uv_your;

std::unique_ptr<ConformalParametrization> C_sys;

// Callback for the GUI
void callback(){
  static int numPoints = 2000;
  static float param = 3.14;

  ImGui::PushItemWidth(100);
  if (ImGui::Button("Show Conformal Parametrization Libigl")) {
    Eigen::VectorXi bnd,b(2,1);
    igl::boundary_loop(F,bnd);
    b(0) = bnd(0);
    b(1) = bnd(bnd.size()/2);
    Eigen::MatrixXd bc(2,2);
    bc<<0,0,1,0;
    //use the in-build Libigl parametrization to compare  
    igl::lscm(V,F,b,bc,V_uv_libigl);


    // Scale the uv
    // V_uv_libigl *= 5;s
    polyscope::getSurfaceMesh("Input Mesh")->addVertexParameterizationQuantity("Conformal Parametrization Libigl", V_uv_libigl);
    polyscope::registerSurfaceMesh2D("Conformal Parametrization Libigl", V_uv_libigl,F);

  }

  if(ImGui::Button("Compute Spectral Conformal Parameterization")){
	C_sys->build_parametrizations();
	V_uv_your = C_sys->V_uv;
	polyscope::registerSurfaceMesh2D("Spectral Conformal Parametrization", V_uv_your,F);
  }

  if (ImGui::Button("Your Conformal Parametrization")) {
    polyscope::getSurfaceMesh("Input Mesh")->addVertexParameterizationQuantity("Your Conformal Parametrization", V_uv_your);
    polyscope::registerSurfaceMesh2D("Spectral Conformal Parametrization", V_uv_your,F);

  }

}

int main(int argc, char *argv[])
{
	std::cout<<" -- Spectral Conformal Parametrization -- "<<std::endl;


	// --- Load the meshes ---
	//------------------------
	if (argc < 2) {
		#ifdef USING_GNU
			std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
			igl::readOFF("../data/camelhead.off", V, F);
			std::string filename = "../data/camelhead.off";
			std::cout << "loading: " << filename << std::endl;
			// Load GNU-specific meshes
		#elif defined(USING_MSVC)
			std::cout << "Using MSVC compiler. Loading MSVC-specific meshes." << std::endl;
			igl::readOFF("../../../../data/camelhead.off", V, F);
			std::string filename = "../../../../../data/camelhead.off";
			std::cout << "loading: " << filename << std::endl;
			// Load MSVC-specific meshes
		#else
			std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
			igl::readOBJ("../data/camelhead.off", V, F);
			std::string filename = "../data/camelhead.off";
			std::cout << "loading: " << filename << std::endl;
			// Load default meshes
		#endif
	}
	else {
		std::string filename = std::string(argv[1]);
		// std::string filename = "../../data/hexagon.off";


		if (filename.size() >= 3 && filename.substr(filename.size() - 3) == "off") {
			std::cout << "loading: " << filename << std::endl;
			// Read the mesh as an OFF file
			igl::readOFF(filename, V, F);
		}

		// Check if the last three letters of the filename are "obj"
		else if (filename.size() >= 3 && filename.substr(filename.size() - 3) == "obj") {
			std::cout << "loading: " << filename << std::endl;
			// Read the mesh as an OBJ file
			igl::readOBJ(filename, V, F);
		}

		// If neither "off" nor "obj", print an error message and exit
		else {
			std::cout << "Error: Unsupported file format. Only .off and .obj formats are supported. " << std::endl;
		}
	}


  

	// //Now set up your own system



	C_sys = std::make_unique<ConformalParametrization>(V,F);

	


	// V_uv_your*=5;

	// --- Visualize the mesh using Polyscope ---
	//--------------------------------------------

	// Options
	polyscope::options::autocenterStructures = true;
	polyscope::view::windowWidth = 1024;
	polyscope::view::windowHeight = 1024;

	// Initialize polyscope
	polyscope::init();

	polyscope::registerSurfaceMesh("Input Mesh", V, F);
  	

	polyscope::state::userCallback = callback;


	// Show the gui
	polyscope::show();

	return 0;
}
