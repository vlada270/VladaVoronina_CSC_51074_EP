// #include <igl/opengl/glfw/Viewer.h>

// std includes
#include <iostream>
#include <ostream>
#include <memory>

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

//Our includes
#include "HalfedgeBuilder.h"
#include "LoopSubdivision.h"
#include "SphereGeneration.h"

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V,V1, V_store;
MatrixXi F,F1, F_store;


void perform_loop_subdivision(MatrixXd &V, MatrixXi &F)
{
	std::unique_ptr<HalfedgeBuilder> builder(new HalfedgeBuilder());
	HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F));  // create the half-edge representation
	std::unique_ptr<LoopSubdivision> loop(new LoopSubdivision(V, F, he));   //
	loop->subdivide();										 // perform one round subdivision
	loop->print(0);

	V = loop->getVertexCoordinates();
	F = loop->getFaces();
	

}

void perform_sphere_generation(MatrixXd &V, MatrixXi &F)
{
	
	std::unique_ptr<HalfedgeBuilder> builder(new HalfedgeBuilder());
	HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
	std::unique_ptr<SphereGeneration> generator(new SphereGeneration(V, F, he));   //
	generator->subdivide();										 // perform one round subdivision
	generator->print(0);

	// // update the current mesh
	V = generator->getVertexCoordinates(); // update vertex coordinates
	F = generator->getFaces();

	
}

void perform_adaptive_loop_subdivision()
{
	
	std::unique_ptr<HalfedgeBuilder> builder(new HalfedgeBuilder());
	HalfedgeDS he = (builder->createMeshWithFaces(V.rows(), F)); // create the half-edge representation
	std::unique_ptr<LoopSubdivision> loop(new LoopSubdivision(V, F, he));   //
	loop->subdivide_adaptive();										 // perform one round subdivision
	loop->print(0);

	// // update the current mesh
	V = loop->getVertexCoordinates(); // update vertex coordinates
	F = loop->getFaces();
}

void save_mesh_to_file()
{
	igl::writeOFF("../data/output.off", V, F);

	#ifdef USING_GNU
			igl::writeOFF("../data/output.off", V, F);
			// Load GNU-specific meshes
	#elif defined(USING_MSVC)
			igl::writeOFF("../../../../data/output.off", V, F);
			// Load MSVC-specific meshes
	#else
			std::cout << "Unknown compiler. Defaulting mesh storing." << std::endl;
			igl::writeOFF("../data/output.off", V, F);
			// Load default meshes
	#endif

}

// Callback for the GUI
void callback() {
	static int numPoints = 2000;
	static float param = 3.14;

	ImGui::PushItemWidth(100);
	if (ImGui::Button("Perform Subdivision")) {
		perform_loop_subdivision(V,F);
        polyscope::registerSurfaceMesh("Subdivision Surface Mesh", V, F);
	}

	if (ImGui::Button("Perform Sphere Generation")) {
		perform_sphere_generation(V,F);
        polyscope::registerSurfaceMesh("Subdivision Surface Mesh", V, F);
	}
	
	if (ImGui::Button("Perform Adaptive Loop Subdivision")) {
		perform_adaptive_loop_subdivision();
		polyscope::registerSurfaceMesh("Subdivision Surface Mesh", V, F);
	}
	// Save To File
	if (ImGui::Button("Save Mesh To File")) {
		save_mesh_to_file();
	}
	//reset the mesh
	if (ImGui::Button("Reset Mesh")) {
		V = V_store;
		F = F_store;
		polyscope::registerSurfaceMesh("Subdivision Surface Mesh", V, F);
	}

	ImGui::PopItemWidth();
}

/**
 * Create a triangle mesh corresponding to an octagon inscribed in the unit circle
 */
void createOctagon(MatrixXd &Vertices, MatrixXi &Faces)
{
	Vertices = MatrixXd(6, 3);
	Faces = MatrixXi(8, 3);

	Vertices << 0.0, 0.0, 1.0,
		1.000000, 0.000000, 0.000000,
		0.000000, 1.000000, 0.000000,
		-1.000000, 0.000000, 0.000000,
		0.000000, -1.000000, 0.000000,
		0.000000, 0.000000, -1.000000;

	Faces << 0, 1, 2,
		0, 2, 3,
		0, 3, 4,
		0, 4, 1,
		5, 2, 1,
		5, 3, 2,
		5, 4, 3,
		5, 1, 4;
}

// ------------ main program ----------------
int main(int argc, char *argv[])
{
	std::cout<<" -- Subdivision Surfaces -- "<<std::endl;


	// --- Load the meshes ---
	//------------------------
	
	// std::string filename = "../../data/letter_a.off";

	if (argc<2) {

		// std::string filename = "../../data/bunny_new.off";
		#ifdef USING_GNU
				std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
				igl::readOFF("../data/letter_a.off", V, F);
				std::string filename = "../data/letter_a.off";
				// Load GNU-specific meshes
		#elif defined(USING_MSVC)
				igl::readOFF("../../../../data/letter_a.off", V, F);
				std::string filename = "../../../../data/letter_a.off";
				// Load MSVC-specific meshes
		#else
				std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
				igl::readOBJ("../data/letter_a.off", V, F);
				std::string filename = "../data/letter_a.off";
				// Load default meshes
		#endif
		
	}
	else {
		std::string filename = std::string(argv[1]);
		// Check if the last three letters of the filename are "off"
		if (filename.size() >= 3 && filename.substr(filename.size() - 3) == "off") {
			// Read the mesh as an OFF file
			igl::readOFF(filename, V, F);
		}
		// Check if the last three letters of the filename are "obj"
		else if (filename.size() >= 3 && filename.substr(filename.size() - 3) == "obj") {
			// Read the mesh as an OBJ file
			igl::readOBJ(filename, V, F);
		}
		// If neither "off" nor "obj", print an error message and exit
		else {
			std::cerr << "Error: Unsupported file format. Only .off and .obj formats are supported." << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}

	V_store = V;
	F_store = F;
	





	// --- Visualize the mesh using Polyscope ---
	//--------------------------------------------

	// Options
	polyscope::options::autocenterStructures = true;
	polyscope::view::windowWidth = 1024;
	polyscope::view::windowHeight = 1024;

	// Initialize polyscope
	polyscope::init();

	polyscope::registerSurfaceMesh("Subdivision Surface Mesh", V, F);

	polyscope::state::userCallback = callback;


	// Show the gui
	polyscope::show();

	return 0;


}
