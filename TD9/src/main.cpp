
//includes LibIGL
#include <iostream>
#include "TrivConnMesh.h"

#include "igl/readOFF.h"
#include "igl/readOBJ.h"

//includes Polyscope
#include "polyscope/polyscope.h"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"

#include <iostream>
#include <unordered_set>
#include <utility>

#include "../../deps/polyscope/deps/args/args/args.hxx"
#include "../../deps/polyscope/deps/json/include/nlohmann/json.hpp"

// The mesh, Eigen representation
Eigen::MatrixXd V;
Eigen::MatrixXi F;

std::unique_ptr<TrivConnMesh> tcods_mesh_ptr;



void callback() {
    if(ImGui::Button("Show Spanning tree")) {
        // visualize the spanning tree for debugging
        polyscope::registerCurveNetwork("spanning tree", tcods_mesh_ptr->V_dual,tcods_mesh_ptr->edges_spanning_tree);
        polyscope::getCurveNetwork("spanning tree")->setEnabled(false);
    }
    if(ImGui::Button("Compute Trivial Connections")) {
        //visualize the parallel transported vector field
        polyscope::getSurfaceMesh("input mesh tcods")->addFaceTangentVectorQuantity("trivial transported vectors", tcods_mesh_ptr->parallel_transported_field_vec,tcods_mesh_ptr->basisX,tcods_mesh_ptr->basisY);
    }
}

int main(int argc, char** argv) {
    // Configure the argument parser
    args::ArgumentParser parser("Trivial Connections On Discrete Surfaces");

    args::Positional<std::string> inFile(parser, "mesh", "input mesh");

    // Options
    polyscope::options::autocenterStructures = true;
    polyscope::view::windowWidth = 1024;
    polyscope::view::windowHeight = 1024;



    // Initialize polyscope
    polyscope::init();

    

    if (argc < 2) {
		#ifdef USING_GNU
				std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
				igl::readOBJ("../data/spot.obj", V, F);
				std::string filename = "../data/spot.obj";
				std::cout << "loading: " << filename << std::endl;
				// Load GNU-specific meshes
		#elif defined(USING_MSVC)
				std::cout << "Using MSVC compiler. Loading MSVC-specific meshes." << std::endl;
				igl::readOBJ("../../../../data/spot.obj", V, F);
				std::string filename = "../../../../../data/spot.obj";
				std::cout << "loading: " << filename << std::endl;
				// Load MSVC-specific meshes
		#else
				std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
				igl::readOBJ("../data/spot.obj", V, F);
				std::string filename = "../data/spot.obj";
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
  
  

  


    // Register the mesh with Polyscope
  

    tcods_mesh_ptr = std::make_unique<TrivConnMesh>(V,F);

    //load the surface mesh
    polyscope::registerSurfaceMesh("input mesh tcods", V, F);

 
  
  
    // Add the callback
    polyscope::state::userCallback = callback;

    // Show the gui
    polyscope::show();

    return 0;
}
