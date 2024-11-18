# TD 5 - Iterative closest points


The provided code skeleton operates in four modes (specified by input arguments 1, 2, 3, or 4) to work with different input meshes. Please ensure that you update the data path in the code according to your system configuration:

    Visual Studio Compilation: The data file is located at "../../../../data".
    Linux, Mac, or Windows Terminal: The data file is located at "../data" (this is the default value in the provided code).

![Interface](../imgs/TD5-interface.png)
The Interface: once each exercise is implemented, it will be accessible via the corresponding botton on the left.


## 1. Iterative closest point (with rotation matrices)

### Some useful functions

Eigen provides the function:

- `m.colwise()` and `m.rowwise()` to apply partial reduction operations to matrices. For example `m.colwise().mean()` returns the column-wise mean of a matrix.
- `v.minCoeff(&min)` returns the maximum of all coefficients and puts in `min` its location.
- `igl::octree` and `igl::knn` compute nearest neighbours in point clouds. Find out more here: [libigl tutorial](https://libigl.github.io/tutorial/) and [knn.h](https://github.com/libigl/libigl/blob/master/include/igl/knn.h)

### ICP with Matrices and SVD

Given two 3d point clouds V1 and V2, we want to compute the rigid transformation that maps V1 to V2, by implementing the *Iterative closest point method with rotation matrices* seen in class. A more detailed algorithm can be found here: [Eggert et al.](http://graphics.stanford.edu/~smr/ICP/comparison/eggert_comparison_mva97.pdf)

#### 1.1 ICP: compute the rigid transformation for corresponded point clouds

- Complete the function `transform` in `ICP.cpp` that modifies V1 to align it to V2 when V1 and V2 are in correspondence.

![Ex1](../imgs/TD5-exo1.png)
**Testing your code**: In order to test your implementation, you can run the program and press the corresponding botton `ICP Step Ex1`.  
See [Eigen Jacobi SVD documentation](https://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html).

#### 1.2 ICP for arbitrary point clouds (no correspondence given)

Now we want to solve the problem above (computing a rigid transformation), for two point clouds whose points are not matched. Feel free to use `igl::knn` and `Eigen/SVD` for this exercise.

- Complete the function `nearest_neighbour` to find, for each point in V1, the nearest neighbour in V2.
- **Convergence**: compute and print (at each iteration) the sum of distances between closest neighbors and check that such distances decrease.
![Ex2](../imgs/TD5-exo2.png)
**Testing your code**: Run the program and press `ICP Ex2` to see an animation of your algorithm. You can also test your implementation on the gargoyle_tri, and egea shapes.

**Bonus**: Implement brute force nearest neighbor and compare it with the octree-based method.

#### 1.3 ICP with Point to Plane

As seen in the lecture, point matching can be difficult if the quality of your point cloud is not good. In such cases, it may be better to minimize the distance from point to the tangent plane.

![Ex3](../imgs/TD5-exo3.png)

Implement the method `nearest_neighbour_point_to_plane` and test your code with `ICP Point to Plane Ex3` botton.

## 2. Normals estimation

### 2.1 Computing normals of point clouds

Given a 3d point cloud V1, we want to compute an estimation of normals as explained in the lecture (using the co-variance matrix).

- Complete the function `k_nearest_neighbour` in `pca.cpp` to find the k nearest neighbours for each point in V1.
- Complete the function `compute_normals` to compute the normals at each point by building a co-variance matrix.

Feel free to use `igl::knn` and `EigenSolver` for this exercise.

**Testing your code**: Run the program using mode 4 which loads points on a sphere as data, and press `PCA Ex4` botton to test your implementation. 
![Normals on a sphere point cloud](../imgs/TD5-exo4.png)
Make sure to enable the visualization of the normals in the interface menu, as in this figure. You are now done with this TD!


