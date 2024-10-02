
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
#include "polyscope/point_cloud.h"

#include "LaplacianMesh.h"



Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd N_faces; 
Eigen::MatrixXd N_vertices; 

Eigen::VectorXd Color;
Eigen::VectorXd heat = Eigen::VectorXd::Zero(V.rows());

std::unique_ptr<LaplacianMesh> lapmeshptr;
std::unique_ptr<Mesh> meshptr;

bool is_animating = false;
bool heatflow = false;
bool use_hedge = false;
bool use_implicit_euler = false;

double tt = 0.0;
float timestep = 0.001;


void updateColor(Eigen::VectorXd &C, double t = 0.0){
	C = Eigen::VectorXd::Zero(V.rows());	
	// TODO: Implement the time-dependent color update logic here
	
}

void updateColorGaussianCurvature(Eigen::VectorXd &C){
	C = meshptr->compute_gaussian_curvature();
}

void callback() {
	static int numPoints = 2000;
	static float param = 3.14;

	ImGui::PushItemWidth(100);

	ImGui::Checkbox("Use HalfEdge", &use_hedge);

	if (ImGui::Button("Update Color")) {
		tt = 0.0;	
		if(is_animating==false){
			is_animating = true;
		}
		else{
			is_animating = false;
		}
	}

	if (ImGui::Button("Compute Mesh statistics")) {
      meshptr->vertexDegreeStatistics();
	}

	if(ImGui::Button("Compute Normals")){
		if(!use_hedge){
			N_faces = meshptr->compute_face_normals();
			N_vertices = meshptr->compute_vertex_normals();
		}else{
			N_faces = meshptr->compute_face_normals_hed();
			N_vertices = -meshptr->compute_vertex_normals_hed();
		}

		polyscope::getSurfaceMesh("Mesh")->addFaceVectorQuantity("Face Normals", N_faces);
		polyscope::getSurfaceMesh("Mesh")->addVertexVectorQuantity("Vertex Normals", N_vertices);
	}

	if (ImGui::Button("Calculate Boundaries")) {
		meshptr->count_boundaries();
	}

	if (ImGui::Button("Visualize Gaussian Curvature")) {
		updateColorGaussianCurvature(Color);
		polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("Gaussian Curvature", Color);
	}

	
	ImGui::Checkbox("Use Implicit Euler", &use_implicit_euler);
	ImGui::SliderFloat("Timestep", &timestep, 0.0f, 0.5f);

	if (ImGui::Button("Start Heat flow")) {
		tt = 0.0;
		
		lapmeshptr->setHeat(heat);
		heatflow = true;
		
	}
		

	if (ImGui::Button("Start Wave equation")) {
		
	}


	if(is_animating){
		updateColor(Color,tt);
		polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("Color", Color)->setMapRange({-5.,5.});
		tt+=0.005;
	}

	if(heatflow){
		if(use_implicit_euler){
			lapmeshptr->heat_step_implicit(heat, double(timestep));
		}else{
			lapmeshptr->heat_step_explicit(heat, double(timestep));
		}
		polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("Heat", heat);
	}

	ImGui::PopItemWidth();

}


	
// ------------ main program ----------------
int main(int argc, char *argv[]) {


	std::cout<<" -- Mesh Datastructure and Laplacian -- "<<std::endl;

	// --- Load the meshes ---
	//------------------------

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
	//initialize the mesh datastructures

	lapmeshptr = std::make_unique<LaplacianMesh>(V, F);
	meshptr = std::make_unique<Mesh>(V, F);

	heat = Eigen::VectorXd::Zero(V.rows());

	// --- Visualize the mesh using Polyscope ---
	//--------------------------------------------

	// Options
	polyscope::options::autocenterStructures = true;
	polyscope::view::windowWidth = 1024;
	polyscope::view::windowHeight = 1024;

	// Initialize polyscope
	polyscope::init();

	polyscope::registerSurfaceMesh("Mesh", V, F);
	

	polyscope::state::userCallback = callback;


	// Show the gui
	polyscope::show();

	return 0;
}
