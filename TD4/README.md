# Meshes
## TD 4 - A first 3D interface to manipulate mesh structures

### 1. Pointer Handling in C++




In this TP on 3D meshes, you will deepen your understanding of pointer management in C++. You’ll work with pointers referencing objects in a dynamic environment, where objects can be modified or deleted. Properly managing pointer updates in such scenarios will be essential for your implementation. By the end of this TD, if you approach it with patience and care, you should feel confident in handling these aspects.




### 2. Geometric Modeling and Data Structures
From a geometric modeling perspective, we will explore two different data structures for triangle mesh manipulation:




- The Halfedge Structure:  
  This structure enables efficient local manipulation of a mesh, such as iterating around a vertex. It offers control over connectivity at the local level.




- Adjacency Matrix Representation:  
  This structure, used in libraries like LibIGL, provides a more global encoding of the mesh using adjacency matrices via vertex labels. While it makes local operations more challenging, such representation offers efficient global computations, which are useful in many applications.




Depending on the context, each structure has its advantages. As you work on this implementation, try to understand the underlying logic of these two connectivity representations. Although specific implementations may vary, the core principles remain consistent across different contexts.




---




## Important Notes




This TD is quite long, but each task is straightforward if you follow the instructions and review any missing C++ concepts. To ensure you have sufficient time, you will have until the day before TD6 to submit your work. Note that this is a mandatory TD. Its result will be used in TD6, TD8 and TD9, but it is independent of TD5 and TD7.








### 1. Data Structures
This TD has two main parts.  




In the first part, you will set up and implement a mesh data structure, including a half-edge structure. The `namespace MeshParts` contains the classes `Vertex`, `PrimalFace`, and `HalfEdge`. In the mesh class header, you’ll find `std::vector`s with shared pointers to vertices, half-edges, and faces.  Handling these pointers requires care, so we provide a C++ pointers overview in the `CPP-Recap` folder.




Regarding the constructor of the mesh, in the method initialize_complex you should create the shared pointers to all the objects that you will be using. Thus, you should create the shared pointers to all the vertices, primal faces and halfedges. You will see that the classes of the namespace MeshParts contain as class attributes itself std::weak_ptr to other mesh parts.
   
In the method compute_half_edges you should now fill the half edge data-structure with life, meaning that you should set all the pointers properly for the mesh parts. We provide you with the following packages to build the HalfEdge data-structure,  `HalfedgeBuilder.cpp` and `HalfedgeDS.cpp` (you don't need to modify them).  Read the header files and cpp files of the `HalfedgeDS.cpp` class to understand how this is working. The general workflow of what you will need to do in this method will look as follows: 
```cpp
void Mesh::compute_half_edges()
{

    //1. Setup Half-Edge Builder
    //A unique pointer to a `HalfedgeBuilder` instance (`builder`) initializes a `HalfedgeDS` data structure. The builder creates a half-edge representation of the mesh using face (`F`) data, facilitating efficient edge operations.
    std::unique_ptr<HalfedgeBuilder> builder(new HalfedgeBuilder());
    HalfedgeDS he = builder->createMesh(V.rows(), F);

    // You can use the provided code to build a half-edge datastructure, that will be purely index based. This means each half edge is represented by an integer. If you call for instance he.getNext(12) you will get the index of the next half edge of the 12th half edge


    //2. Half-Edge Creation and Initialization: 
    // After construting HalfedgeDS he, you should create according to the indices of he instances of MeshParts::HalfEdge and store shared pointers to these instances in the std::vector<HalfEdgePtr> this->hedges.

    //3. Half-Edge Connection:
    //Now, that the MeshParts::HalfEdge instances are created, you can use the methods of HalfedgeDS ( read the header file to see what methods are available) to store on each half edge the shared pointers to the start vertex, end vertex, next half edge, flip half edge, boundary etc.

    //4. Vertex and Face Linking
    //  Also, each instance of PrimalFace and Vertex that you created in the method  initialize complex has as attributes std::vectors to pointers of halfedges. Now that the half edges are initialized
    

    //5.  Dual Face Creation (Optional):
   //This segment, partly commented out, would generate dual edges if activated, representing edges that connect the circumcenters of neighboring faces in the dual mesh. It’s useful in advanced mesh structures where dual topology is required. 

   // Let's go! 
}

```
We highly recommend using the visual interface of polyscope to debug your construction, by clicking on different faces for instance and visualizing the corresponding indices, while checking their consistency. The methods  of HalfedgeDS will primarily return indices to other halfedges, faces, vertices.




A rich set of tests is also provided for this TD, that you can use to test your implementation once you have all the data-structure nicely defined. The first part of the unittests will basically check whether you set all the pointers properly.

**Remark:**
Once you try to calculate the mesh datastructure for an open mesh, i.e a mesh with boundary, you need to be careful when you set the pointers for HalfEdges on the boundary. For your own sanity, you could add a condition in the method `HalfEdge::getFlipHalfEdge()` that will check (and throw a warning) if you try to call this method for a boundary half edge. In general, feel free to add more attributes or methods to the mesh parts if it helps you in your implementation.   

### 1.1 Vertex statistics
Once you are done with the constructor of the class `Mesh`, to get familiar with the new datastructure compute the degree of all vertices at once using the face-vertex structure of LibiGL. Complete the method 
`vertexDegreeStatistics`.
Use your HalfEdge datastructure as a test for yourself to see if the degrees are consistent.  
You can use `std::chrono::high_resolution_clock::now()` to measure performance. Measure the performance if you only calculate the degree for a single vertex, rather then all of the vertices. 

### 1.2 Time Dependent Coloring
To change the color of the mesh over time, update the color in each frame. In `main.cpp` complete the method
```cpp
void updateColor(Eigen::VectorXd &C, double t = 0.0){
    C = Eigen::VectorXd::Zero(V.rows());    
    // TODO: Implement the time-dependent color update logic here
   
}
```
where you should implement a time dependant color function. Once you tick the checkbox `is animating` a boolean variable is set and you will find in the callback method
```cpp
if(is_animating){
    updateColor(Color,tt);
    polyscope::getSurfaceMesh("Mesh")->addVertexScalarQuantity("Color", Color)->setMapRange({-5.,5.});
    tt+=0.005;
}
```
You do not need to edit this method, but for your understanding what happens here is the following:




1) We load a mesh in the `main.cpp` with for example `igl::readOBJ` that converts into a `MatrixXd V` and `MatrixXi F` that contains the point positions and mesh information.




2) Next, we want to be able to visualize the mesh with Polyscope. In Polyscope you need to register an instance of a surface structure, this can be done with
```cpp
 polyscope::registerSurfaceMesh("Mesh", V, F);
```
Together with `polyscope::show();` this allows to show, with the help of the openGL backend, the mesh. Further, the method `polyscope::state::userCallback = callback;` will basically in each frame check if there is some change on the mesh or related information and walk through the `callback` function. This is where it will find the aforementioned lines for the update of the color  

In your tests experiment with the function `update_color`.


### 1.3 Normal Computation
The goal is to compute the normals using different methods. In the class `Mesh` you will find the methods
```cpp
Eigen::MatrixXd compute_vertex_normals();
Eigen::MatrixXd compute_vertex_normals_hed();

Eigen::MatrixXd compute_face_normals();
Eigen::MatrixXd compute_face_normals_hed();
```
where you should calculate the vertex and face normals using the face-based datastructure and the half-edge datastructure. You will see in the `struct` for `Vertex` and `PrimalFace` that there is an attribute for the normal. Overwrite this attribute when calling the methods to compute normals. Once they are implemented you can check with the gui and the button `Compute Normals` whether your normal field looks good. You can visualize the results of both implementations via the interface options as in this figure.

![Normal Computation](../imgs/normals_bunny.png)

### 1.4 Number of Boundaries
Compute the number of boundary components in the mesh (complete `countBoundaries`).

### 1.5 Gaussian Curvature Computation
We will use the half edge data structure to compute the Gaußian curvature at each vertex. Complete the function `compute_gaussian_curvature`. Normalize based on the area of the dual (Voronoi) cell. You can visualize the Gaußian curvature as a scalar-valued function by clicking on the corresponding button.

![Gaußian curvature on the torus](../imgs/torus_curvature.png)

## 2. Mesh Operators


As you have seen in class, the Laplace-Beltrami operator on a mesh is defined as a matrix $\Delta = A^{-1}L \in \mathbb{R}^{\mathcal{V}\times\mathcal{V}}$, where




$$ L_{ij} = \begin{cases} \frac{1}{2}( \text{cot}(\alpha_{ij})+\text{cot}(\beta_{ij})) & \text{if } i \in \mathcal{N}(j)\\ -\sum_{k \in \mathcal{N}(i),\ k\neq i} L_{ik} & \text{if } i = j \\ 0 & \text{if } i \notin \mathcal{N}(j) \end{cases}$$
and $A = \mathrm{diag}((\mathrm{Voronoi-Area})_i)$.
We will now in the following build on top of the mesh datastructure that you implemented before also construct the Laplace-Beltrami operator. We construct therefore a novel class `LaplacianMesh` that inherits from the class `Mesh`.
#### Sparse matrices in Eigen
You will notice that for an ordinary mesh, almost all entries of the Laplace-matrix are zero. Therefore, instead of storing a huge matrix that contains mainly zeros, we will instead use a datastructure, where we only store the non-zero entries. These matrices are called Sparse Matrices. You will find in the file `LaplacianMesh.h` a `typedef Eigen::SparseMatrix<double> SpMat;`, meaning that the Eigen Sparse matrices whose entries are `double` can be identified with this call. A sparse matrix with double values can be declared through Eigen::SparseMatrix<double>. In order to create such a matrix we need to create Triplets
`Eigen::Triplet<double>(row, column, value)` and stack them in a vector. Then we can use the method `setFromTriplets` to fill the sparse matrix. For these kind of matrices, Eigen has special efficient in-build solvers that we are going to use in the following.




Complete the methods needed for the constructor of the class `LaplacianMesh`.


To verify your implementation, you will will now implement a UnitTest comparing your results against the pre-build methods from Libigl.


#### Unit Tests with the Googletest Library - Step by Step Guide

In previous exercises, we've introduced the concept of unit testing. For this, we use the Googletest library. The basic idea is to create a separate executable that runs a series of tests, checking whether the actual outcomes match expected results defined at the start. In the `CMakeLists.txt`, you’ll find the following lines:

```cpp
message("\n\n == CMAKE including the googletest\n")
add_subdirectory(../deps/googletest deps/googletest)

add_executable(unit_tests tests/unittests.cpp)
target_link_libraries(unit_tests td_4_lib gtest gtest_main)
```

This configuration links the test executable to both our custom library (e.g., mesh files) and the Googletest library. However, note that other libraries, like Polyscope, may not be needed for verifying certain functionalities for mesh correctness. While it’s possible to include the Polyscope header without errors during compilation, issues will arise during the linking phase, causing the build process to fail. This is because the linker won’t be able to resolve the necessary references.

You will now check, whether your implementation of the matrix `L` matches the pre-computed cotangent matrix of Libigl. Therefore, go to the file `unittests.cpp` and include the header `#include "igl/cotmatrix.h"` Also, include the header for your class `LaplacianMesh.h` in the `unittests.cpp`. Next, inside the test file create an instance of your class `LaplacianMesh`. For example 
`LaplacianMesh LapMesh= LaplacianMesh(V0,F0);`

Next, inside the Test for the cotangent matrix call the Libigl method for the construction of the cotangent matrix. For example

```cpp
TEST(LaplacianTest,CheckWithLibigl) {
  SpMat L_check;
  igl::cotmatrix(V0,F0,L_check);
}
```

You have now created the test values for your matrix L. Next, we will iterate over the non-zero entries of matrix L and verify that they match the precomputed values. The goal is to iterate over the non-zero elements of your sparse matrix and compare them with the corresponding entries in the libigl matrix. This can be done as follows:

```cpp
 for (int k = 0; k < L_check.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(L_check, k); it; ++it) {
          int i = it.row();
          int j = it.col();
          double entry_libigl_ij = it.value();
          double your_entry_ij = LapMesh.L.coeffRef(i,j);
      }
    }

```

Now, that you have your value and the comparison value, you can indeed do the googletest assertion. There are multiple checks that can be done. A full overview can be found [here](https://google.github.io/googletest/reference/assertions.html). For the test here, we want to check if two double values match. Here, for instance the assertion `EXPECT_NEAR` is what you want. The full example could look as follows

```cpp
#include "../src/LaplacianMesh.h"
LaplacianMesh LapMesh = LaplacianMesh(V0,F0);
TEST(LaplacianTest,CheckWithLibigl) {
  SpMat L_check;
  igl::cotmatrix(V0,F0,L_check);
  for (int k = 0; k < L_check.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(L_check, k); it; ++it) {
      int i = it.row();
      int j = it.col();
      double entry_libigl_ij = it.value();
      double your_entry_ij = LapMesh.L.coeffRef(i,j);
      EXPECT_NEAR(entry_libigl_ij,your_entry_ij,1e-8);
      
    }
  }
}

```

Now, repeat the same for the second test with the area matrix, you can use `igl::massmatrix(V0, F0, igl::MASSMATRIX_TYPE_VORONOI, A_check);` to compare.



### 2.1 The Heat Flow
The heat equation is given by:


$$ \frac{\partial u}{\partial t} = \Delta u $$


Given the heat at time $u^{k}$, you can calculate heat at time $u^{k+1}$ using:




Explicit integration: $u^{k+1} = dt \cdot A^{-1} \cdot L \cdot u^{k} + u^{k}$


Semi-implicit integration: $(A - dt \ L) u^{k+1} = A u^{k}$


Implement both in `heat_step_explicit` and `heat_step_implicit`. Note that for the implicit Euler, you will need to solve in each step a Linear System. Fortunately Eigen provides us with a set of tailored sparse solvers.
The steps to use these solvers is as follows, illustrated in the example of a conjugate gradient solver. 

First, create an instance of a Sparse solver
```cpp
Eigen::ConjugateGradient<SpMat, Eigen::Upper> solver;
```
To solve the system $A u = b$, where $A$ is a sparse matrix, you can now pass the sparse matrix to the solver and solve through 
```cpp
solver.compute(A);
u = solver.solve(b);
```
 An overview can be found [here](https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html) and depending on the structure of your linear system, different solvers might be appropriate. 

 Visualize heat flow by modifying the animation loop. Experiment how big you can choose the step size of your update scheme to get a stable solution.

 <img src="../imgs/heat_bunny.png" width="500" height="500" />



**Bonus:**

The wave equation is given by
$$\frac{\partial^2 u}{\partial t^2} = \Delta u$$


Set appropriate initial conditions for the elongation and velocity, update the time integrator for a second order scheme and visualize surface waves.





