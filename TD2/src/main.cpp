// std includes
#include <iostream>
#include <ostream>
#include <memory>
#include <chrono>
#include <vector>
#include <thread>
#include <string>

// libigl includes
#include "igl/readPLY.h"


// Eigen includes
#include "Eigen/Core"

//Polyscope includes
#include "polyscope/polyscope.h"
#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

#include "util.h"
#include "lagrange.cpp"
#include "linear.cpp"
#include "cubic.cpp"
#include "hermite.cpp"

using Eigen::MatrixXi, Eigen::MatrixXd, Eigen::Vector3d;
MatrixXd V1; // matrix storing vertex coordinates of the input curve
MatrixXi F1;
bool linear = true;
bool cubic = false;
bool lagrange = false;
bool hermite = false;

bool cubic_tangents = false;
bool hermite_tangents = false;

int resolution = 5; // resolution of displayed curve

LinearInterpolation interp_linear;
LagrangeInterpolation interp_lagrange;
CubicInterpolation interp_cubic;
HermiteInterpolation interp_hermite;


polyscope::PointCloud* psPoints;



void callback(){
  static int numPoints = 2000;
	static float param = 3.14;

	ImGui::PushItemWidth(100);

  if(ImGui::SliderInt("Resolution", &resolution, 2, 2000)){}

  if (ImGui::Checkbox("Linear", &linear)) {
    cubic = false;lagrange = false;hermite = false;linear = true;
  }
  if (ImGui::Checkbox("Lagrange", &lagrange)) {
    linear = false;cubic = false;hermite = false;lagrange = true;
  }

  if (ImGui::Checkbox("Cubic", &cubic)) {
    linear = false;lagrange = false;hermite = false;cubic = true;
  }
  ImGui::SameLine();
  if(ImGui::Checkbox("Cubic Tangents", &cubic_tangents)){}

  if (ImGui::Checkbox("Hermite", &hermite)) {
    linear = false;cubic = false;lagrange = false;hermite = true;
  }
  ImGui::SameLine();
  if(ImGui::Checkbox("Hermite Tangents", &hermite_tangents)){}

  if(ImGui::Button("Interpolate")){

    if(linear){
       
      MatrixXd linspace = MatrixXd::Zero(resolution, 3);
      util::build_linspace(linspace, V1); // initialize the X axis of the interpolation
      std::vector<Vector3d> points;
      std::vector<std::array<int, 2>> edges;
      std::vector<std::array<double, 3>> colors;

      for (size_t i = 0; i < resolution; i++) {
        double time = linspace(i, 0);
        linspace(i, 1) = interp_linear.eval_function(time);
        Vector3d pt = Vector3d(linspace(i, 0), linspace(i, 1), 0);
        points.push_back(pt);
        if(i<resolution-1){
          edges.push_back({int(i), int(i+1)});
          colors.push_back({0.0, 0.0, 1.0});
        }
      }

      polyscope::registerCurveNetwork("Linear Interpolation", points, edges);

      

    }

    if(lagrange){
      MatrixXd linspace = MatrixXd::Zero(resolution, 3);
      util::build_linspace(linspace, V1); // initialize the X axis of the interpolation
      std::vector<Vector3d> points;
      std::vector<std::array<int, 2>> edges;
      std::vector<std::array<double, 3>> colors;

      for (size_t i = 0; i < resolution; i++) {
        double time = linspace(i, 0);
        linspace(i, 1) = interp_lagrange.eval_function(time);
        Vector3d pt = Vector3d(linspace(i, 0), linspace(i, 1), 0);
        points.push_back(pt);
        if(i<resolution-1){
          edges.push_back({int(i), int(i+1)});
          colors.push_back({0.0, 0.0, 1.0});
        }
      }

      polyscope::registerCurveNetwork("Lagrange Interpolation", points, edges);
      polyscope::getCurveNetwork("Lagrange Interpolation")->resetTransform();
      polyscope::getPointCloud("Input Points")->resetTransform();
    }

    if(cubic){
      MatrixXd linspace = MatrixXd::Zero(resolution, 3);
      util::build_linspace(linspace, V1); // initialize the X axis of the interpolation
      std::vector<Vector3d> points;
      std::vector<std::array<int, 2>> edges;
      std::vector<std::array<double, 3>> colors;

      for (size_t i = 0; i < resolution; i++) {
        double time = linspace(i, 0);
        linspace(i, 1) = interp_cubic.eval_function(time);
        Vector3d pt = Vector3d(linspace(i, 0), linspace(i, 1), 0);
        points.push_back(pt);
        if(i<resolution-1){
          edges.push_back({int(i), int(i+1)});
          colors.push_back({0.0, 0.0, 1.0});
        }
      }

      polyscope::registerCurveNetwork("Cubic Interpolation", points, edges);

      if(cubic_tangents){
        std::vector<Vector3d> tangent_vectors;
        tangent_vectors.push_back(Vector3d(0, 0, 0));
        for (size_t i = 1; i < V1.rows()-1; i++) {
          MatrixXd dX(V1.rows(), 2);
          interp_cubic.eval_tangent(i, dX, V1(i, 0));
          Vector3d tangent = Vector3d(dX(i, 0), dX(i, 1), 0);
          tangent_vectors.push_back(tangent);
          
        }
        tangent_vectors.push_back(Vector3d(0, 0, 0));

        polyscope::getPointCloud("Input Points")->addVectorQuantity("Tangents", tangent_vectors);
      }
    }


    if(hermite){
      interp_hermite = HermiteInterpolation(V1);
      MatrixXd linspace = MatrixXd::Zero(resolution, 3);
      util::build_linspace(linspace, V1); // initialize the X axis of the interpolation
      std::vector<Vector3d> points;
      std::vector<std::array<int, 2>> edges;
      std::vector<std::array<double, 3>> colors;
      interp_hermite.eval_function(linspace);

      for(size_t i = 0; i < linspace.rows(); i++){
        Vector3d pt = Vector3d(linspace(i, 0), linspace(i, 1), 0);
        points.push_back(pt);
        if(i<resolution-1){
          edges.push_back({int(i), int(i+1)});
          colors.push_back({0.0, 0.0, 1.0});
        }
      }

      polyscope::registerCurveNetwork2D("Hermite Interpolation", points, edges);

      if(hermite_tangents){

        //TODO
      }
    }


  }
}




int main(int argc, char *argv[])
{

    // --- Load the meshes ---
	//------------------------
	

    if(argc<2) {
        

        #ifdef USING_GNU
                std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
                igl::readPLY("../data/curve.ply", V1, F1);
                std::string filename = "../data/curve.ply";
                std::cout << "loading: " << filename << std::endl;
                // Load GNU-specific meshes
        #elif defined(USING_MSVC)
                std::cout << "Using MSVC compiler. Loading MSVC-specific meshes." << std::endl;
                igl::readPLY("../../../../data/curve.ply", V1, F1);
                std::string filename = "../../../../../data/curve.ply";
                std::cout << "loading: " << filename << std::endl;
                // Load MSVC-specific meshes
        #else
                std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
                igl::readPLY("../data/curve.ply", V, F);
                std::string filename = "../data/curve.ply";
                std::cout << "loading: " << filename << std::endl;
                // Load default meshes
        #endif

    }
    else {
        std::string filename = std::string(argv[1]);
        if (filename.size() >= 3 && filename.substr(filename.size() - 3) == "ply") {
            // Read the mesh as an OFF file
            igl::readPLY(filename, V1, F1);
            std::cout << "reading input file: " << filename << std::endl;
}
        else {
            std::cout << "Error: input file must be a .ply file" << std::endl;
            return 0;
        }
    }
  
  
  //  print the number of mesh elements
  std::cout << "Points: " << V1.rows() << std::endl;




  std::cout<<" -- Interpolation of Curves -- "<<std::endl;

  // // choose interpolation method
  interp_linear = LinearInterpolation(V1);
  interp_lagrange = LagrangeInterpolation(V1);
  interp_cubic = CubicInterpolation(V1);
  


	// --- Visualize the mesh using Polyscope ---
	//--------------------------------------------

	// Options
	polyscope::options::autocenterStructures = true;
	polyscope::view::windowWidth = 1024;
	polyscope::view::windowHeight = 1024;

	// Initialize polyscope
	polyscope::init();

	psPoints = polyscope::registerPointCloud("Input Points", V1);
	

	polyscope::state::userCallback = callback;


	// Show the gui
	polyscope::show();

	return 0;



  
  return(0);
}

