// std includes
#include <iostream>
#include <ostream>
#include <memory>
#include <chrono>
#include <vector>
#include <thread>
#include <array>
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

#include "Bezier.h"

using namespace Eigen;

Eigen::MatrixXd V1; // matrix storing vertex coordinates of the input curve
Eigen::MatrixXi F1; // matrix storing face indices of the input curve
bool is_animating = false;
double timer = 0.0;
int levels = 3;

double initial_radius = 0.1;
int subdivisions_ring = 10;
polyscope::CurveNetwork* psCurve;
polyscope::PointCloud* psPoints;



namespace util{

  void build_linspace(MatrixXd &linspace, const MatrixXd &V)
  {
    for (size_t i = 0; i < linspace.rows(); i++)
    {
      linspace(i, 0) = V.col(0).minCoeff() + ((V.col(0).maxCoeff() - V.col(0).minCoeff()) / (linspace.rows() - 1)) * i;
    }
  }



  void add_curve(const MatrixXd &V, std::vector<Vector3d>& points, std::vector<std::array<int, 2>>& edges,std::vector<std::array<double, 3>>& colors, std::array<double, 3> color = {0.0, 0.0, 0.0}){
    // Convert points to a suitable format for Polyscope
    int n_pts = points.size();
    for(int i = 0; i < V.rows(); i++){
      points.push_back(Eigen::Vector3d(V(i, 0), V(i, 1), V(i, 2)));
      if(i<V.rows()-1){
        edges.push_back({n_pts+i, n_pts + (i+1)});
        colors.push_back(color);
      }
    }
  }
  

  
  void translate_points(MatrixXd &V, MatrixXd a){
    for(int i = 0; i < V.rows(); i++){
      V.row(i) += a;
    }
  }
  
}


void test1()
{
  std::cout << "basic rendering of a Bezier curve of degree " << (V1.rows() - 1) << std::endl;
  Bezier bezier;
  int resolution = 100;
  std::cout<<bezier.de_casteljau(V1,0.76358)<<std::endl;
  MatrixXd plot = bezier.plot_curve(V1, resolution);
  std::vector<Vector3d> pts_for_curve_network;
  std::vector<std::array<int, 2>> edges_for_network;
  std::vector<std::array<double, 3>> colors;
  util::add_curve(plot, pts_for_curve_network, edges_for_network,colors);
  util::add_curve(V1, pts_for_curve_network, edges_for_network,colors, {1.0, 0.0, 0.0});
  // util::draw_control_polygon(V1);
  psCurve = polyscope::registerCurveNetwork("Bezier Curve", pts_for_curve_network,edges_for_network);
  polyscope::getCurveNetwork("Bezier Curve")->addEdgeColorQuantity("Coloring Curve", colors)->setEnabled(true);
}

void test2(double time_step){

  // std::cout<<"Visualize the intermediate steps of the Algorithm"<<std::endl;

  Bezier bezier;
  
  std::vector<Vector3d> pts_for_curve_network;
  std::vector<std::array<int, 2>> edges_for_network;
  std::vector<std::array<double, 3>> colors;

  util::add_curve(V1, pts_for_curve_network,edges_for_network,colors, {1.0, 0.0, 0.0});
  MatrixXd V2 = V1;
  while(V2.rows()>1){
    V2 = bezier.de_casteljau_intermediate(V2,time_step);
    util::add_curve(V2, pts_for_curve_network,edges_for_network,colors, {0.0, 0.0, 1.0});
  }

  //now draw on the same data channel the curve

  MatrixXd pts_for_plot = bezier.plot_curve(V1,100,time_step);

  util::add_curve(pts_for_plot, pts_for_curve_network,edges_for_network,colors, {0.0, 1.0, 0.0});

  psCurve = polyscope::registerCurveNetwork("Bezier Curve", pts_for_curve_network,edges_for_network);
  polyscope::getCurveNetwork("Bezier Curve")->addEdgeColorQuantity("Coloring Curve", colors)->setEnabled(true);
  psPoints = polyscope::registerPointCloud("Control Points", V1);
}

void test3()
{
  std::cout << "splitting a Bezier curve of degree " << (V1.rows() - 1) << std::endl;
  Bezier bezier;
  int resolution = 500;
  std::array<double,3> green = {0.1, 0.9, 0.1};
  std::array<double,3> yellow = {1.0, 1.0, 0.1};

  std::vector<Vector3d> pts_for_curve_network;
  std::vector<std::array<int, 2>> edges_for_network;
  std::vector<std::array<double, 3>> colors;

  std::vector<MatrixXd> curves = bezier.subdivide(V1, 0.5); // subdivide the curve for t=0.5
  MatrixXd b0 = curves[0];
  MatrixXd b1 = curves[1];

  MatrixXd p(1, 3), q(1, 3);
  p(0, 0) = 1.0; p(0, 1) = 0.0; p(0, 2) = 1.0;
  q(0, 0) = 1.0; q(0, 1) = 0.0; q(0, 2) = -0.8;

  MatrixXd plot = bezier.plot_curve(V1, resolution); //the input curve B(t)

  util::add_curve(plot, pts_for_curve_network,edges_for_network,colors, {0.0, 0.0, 0.0});
  

  // b0 = curve_renderer.translate_points(b0, p);
  util::translate_points(b0, p);

  MatrixXd plot0 = bezier.plot_curve(b0, resolution); // draw the first sub-curve

  util::add_curve(plot0, pts_for_curve_network,edges_for_network,colors, green);
  

  util::translate_points(b1, q);
  MatrixXd plot1 = bezier.plot_curve(b1, resolution); // draw the second sub-curve
  util::add_curve(plot1, pts_for_curve_network,edges_for_network,colors, yellow);
  // curve_renderer.draw_colored_curve(plot1, yellow);
  // //curve_renderer.draw_control_polygon(b1);

  psCurve = polyscope::registerCurveNetwork("Bezier Curve", pts_for_curve_network,edges_for_network);
  polyscope::getCurveNetwork("Bezier Curve")->addEdgeColorQuantity("Coloring Curve", colors)->setEnabled(true);
  
}

void test4()
{
  std::cout << "rendering a Bezier curve with recursive subdivision" << std::endl;
  Bezier bezier;
  
  MatrixXd bounding_box = bezier.subdivision_plot(V1, 0);

  MatrixXd plot = bezier.subdivision_plot(V1, levels);
  
  std::vector<Vector3d> pts_for_curve_network;
  std::vector<std::array<int, 2>> edges_for_network;
  std::vector<std::array<double, 3>> colors;

  util::add_curve(bounding_box, pts_for_curve_network,edges_for_network,colors, {0.0, 0.0, 0.0});
  util::add_curve(plot, pts_for_curve_network,edges_for_network,colors, {1.0, 0.0, 0.5});

  psCurve = polyscope::registerCurveNetwork("Bezier Curve", pts_for_curve_network,edges_for_network);
  polyscope::getCurveNetwork("Bezier Curve")->addEdgeColorQuantity("Coloring Curve", colors)->setEnabled(true);


}

void test5()
{
  std::cout << "Computing tangents and normals " << std::endl;
  Bezier bezier;
  int resolution = 500;

  MatrixXd plot = bezier.plot_curve(V1, resolution);
  
  std::vector<Vector3d> pts_for_curve_network;
  std::vector<std::array<int, 2>> edges_for_network;
  std::vector<std::array<double, 3>> colors;

  util::add_curve(plot, pts_for_curve_network,edges_for_network,colors, {0.0, 0.0, 0.0});

  // util::add_curve(V1, pts_for_curve_network,edges_for_network,colors, {1.0, 0.0, 0.0});
  
  std::vector<Vector3d> evaluation_points;
  std::vector<Vector3d> tangents;
  std::vector<Vector3d> normals;

  for (double step = 0.1; step < 1.; step = step + 0.1)
  {
    Vector3d a = bezier.toVector(bezier.de_casteljau(V1, step));
    evaluation_points.push_back(a);
    MatrixXd tangent = bezier.toVector(bezier.compute_tangent(V1, step));
    tangents.push_back(tangent);
    MatrixXd normal = bezier.toVector(bezier.compute_normal(V1, step));
    normals.push_back(normal);
  }

  psCurve = polyscope::registerCurveNetwork("Bezier Curve", pts_for_curve_network,edges_for_network);
  polyscope::getCurveNetwork("Bezier Curve")->addEdgeColorQuantity("Coloring Curve", colors)->setEnabled(true);
  polyscope::registerPointCloud("Evaluation Points", evaluation_points);
  polyscope::getPointCloud("Evaluation Points")->addVectorQuantity("Tangents", tangents)->setEnabled(true);
  polyscope::getPointCloud("Evaluation Points")->addVectorQuantity("Normals", normals)->setEnabled(true);
}

void test6()
{
  std::cout << "3D rendering of a Bezier curves " << std::endl;
  Bezier bezier;
  int resolution = 500; // for rendering the bezier curve
  int w=8; // width of the mesh
  int h=100; // height of the mesh
  double delta=1.0/h;

  MatrixXd plot = bezier.plot_curve(V1, resolution);

  std::vector<Vector3d> pts_for_curve_network;
  std::vector<std::array<int, 2>> edges_for_network;
  std::vector<std::array<double, 3>> colors;

  util::add_curve(plot, pts_for_curve_network,edges_for_network,colors, {0.0, 0.0, 0.0});
  //curve_renderer.draw_control_polygon(V1);

  
  std::vector<MatrixXd> loops_for_network;
  double t=0.0;
  for (int i=0; i<h; i++)
  {
    MatrixXd loop = bezier.compute_loop_of_vertices(V1, t, subdivisions_ring, 2.0*initial_radius);
    loops_for_network.push_back(loop);
    t=t+delta;
    initial_radius+=delta;
    
    util::add_curve(loop, pts_for_curve_network,edges_for_network,colors, {1.0, 0.0, 0.0});

  }
  psCurve = polyscope::registerCurveNetwork("Bezier Curve", pts_for_curve_network,edges_for_network);
  polyscope::getCurveNetwork("Bezier Curve")->addEdgeColorQuantity("Coloring Curve", colors)->setEnabled(true);
}

void callback(){
  static int numPoints = 2000;
	static float param = 3.14;

	ImGui::PushItemWidth(100);

  

	if (ImGui::Button("Task 1")) {
		test1();
	}


  
  
  if (ImGui::Button("Task 2")) {
    if(!is_animating){
      test2(timer);
    }
    
	}
  ImGui::SameLine();
  if (ImGui::InputDouble("Time Step", &timer)){}
  ImGui::SameLine();
  ImGui::Checkbox("Animating", &is_animating);

  if (ImGui::Button("Task 3")) {
    test3();
  }

  ImGui::SameLine();
  if (ImGui::InputInt("Number of Recursion levels", &levels)) {}

  if(ImGui::Button("Task 4")){
    test4();
  }

  if(ImGui::Button("Task 5")){
    test5();
  }

  if(ImGui::Button("Task 6")){
    test6();
  }
  ImGui::SameLine();
  if(ImGui::InputDouble("Initial Radius", &initial_radius)){}
  ImGui::SameLine();
  if(ImGui::InputInt("Subdivisions Ring", &subdivisions_ring)){}
  
  

  if(is_animating){
    if(timer<1.0){
      timer += 0.05;
      test2(timer);
    }else{
      timer = 0.0;
      is_animating = false;
      test2(timer);
    }
    
    // Sleep for 500 milliseconds
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    
  }

  if(timer==1.0){
    timer = 0.;
    is_animating = false;
    test2(timer);
  }
}

int main(int argc, char *argv[])
{

    // --- Load the meshes ---
    //------------------------
    
    if(argc<2) {
        #ifdef USING_GNU
                std::cout << "Using GNU compiler. Loading GNU-specific meshes." << std::endl;
                igl::readPLY("../data/cubic.ply", V1, F1);
                std::string filename = "../data/cubic.ply";
                std::cout << "loading: " << filename << std::endl;
                // Load GNU-specific meshes
        #elif defined(USING_MSVC)
                std::cout << "Using MSVC compiler. Loading MSVC-specific meshes." << std::endl;
                igl::readPLY("../../../../data/cubic.ply", V1, F1);
                std::string filename = "../../../../../data/cubic.ply";
                std::cout << "loading: " << filename << std::endl;
                // Load MSVC-specific meshes
        #else
                std::cout << "Unknown compiler. Defaulting mesh loading." << std::endl;
                igl::readPLY("../data/cubic.ply", V1, F1);
                std::string filename = "../data/cubic.ply";
                std::cout << "loading: " << filename << std::endl;
                // Load default meshes
        #endif
    }
    else {
        std::string filename = std::string(argv[1]);
        std::cout << "reading input file: " << filename << std::endl;
        igl::readPLY(filename, V1, F1);
    }
    
  
    double time = 0.3;

    //  print the number of mesh elements
    std::cout << "Points: " << V1.rows() << std::endl;

    // igl::opengl::glfw::Viewer viewer; // create the 3d viewer
    // Util curve_renderer(viewer);      // 3D renderer for drawing points, polylines, curves, ...


    std::cout<<" -- Bezier Curves -- "<<std::endl;

	


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


}
