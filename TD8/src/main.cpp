#include <iostream>

#include "igl/readOFF.h"
#include "igl/readOBJ.h"

//includes Polyscope
#include "polyscope/polyscope.h"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"

#include <iostream>
#include <unordered_set>
#include <utility>

#include "../../deps/polyscope/deps/args/args/args.hxx"
#include "../../deps/polyscope/deps/json/include/nlohmann/json.hpp"

#include "ElectricMesh.h"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
MatrixXi F;

std::unique_ptr<ElectricMesh> test_electric_mesh;

// This function is called every time a keyboard button is pressed

void callback() {
  if(ImGui::Button("Solve For Electric Field")) {
    
	MatrixXd rho_in = MatrixXd::Zero(V.rows(), 1);
	
	test_electric_mesh->initialize_charge_density(rho_in);
	test_electric_mesh->solve_for_u();
	test_electric_mesh->compute_electric_field();
	
	//visualize the charge density
    polyscope::getSurfaceMesh("input electric mesh")->addVertexScalarQuantity("charge density", test_electric_mesh->rho_vector);
    //visualize the electric field 
    polyscope::getSurfaceMesh("input electric mesh")->addFaceTangentVectorQuantity("electric field", test_electric_mesh->electric_field_in_frame,test_electric_mesh->basisX,test_electric_mesh->basisY);
  }
}


// ------------ main program ----------------
int main(int argc, char *argv[])
{

	// Configure the argument parser
  args::ArgumentParser parser("Dual Complex And Hodge Helmholtz Decomposition");

  args::Positional<std::string> inFile(parser, "mesh", "input mesh");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Options
  polyscope::options::autocenterStructures = true;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;



  // Initialize polyscope
  polyscope::init();

  

  if (argc < 2) {
		#ifdef USING_GNU
				std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
				igl::readOFF("../data/bunny_new.off", V, F);
				std::string filename = "../data/bunny_new.off";
				std::cout << "loading: " << filename << std::endl;
				// Load GNU-specific meshes
		#elif defined(USING_MSVC)
				std::cout << "Using MSVC compiler. Loading MSVC-specific meshes." << std::endl;
				igl::readOFF("../../../../data/bunny_new.off", V, F);
				std::string filename = "../../../../../data/bunny_new.off";
				std::cout << "loading: " << filename << std::endl;
				// Load MSVC-specific meshes
		#else
				std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
				igl::readOFF("../data/bunny_new.off", V, F);
				std::string filename = "../data/bunny_new.off";
				std::cout << "loading: " << filename << std::endl;
				// Load default meshes
		#endif
	}
	else {
		std::string filename = std::string(argv[1]);
		// std::string filename = "../../data/bunny_new.off";


		// Check if the last three letters of the filename are "off"
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
	}
  
  


	test_electric_mesh = std::make_unique<ElectricMesh>(V,F);
	

	polyscope::init(); // Initialize polyscope
	//load the surface mesh
	polyscope::registerSurfaceMesh("input electric mesh", test_electric_mesh->V, test_electric_mesh->F);



	// Add the callback
	polyscope::state::userCallback = callback;

	// Show the gui
	polyscope::show();

	return 0;

	
}
