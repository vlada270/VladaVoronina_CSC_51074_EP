#ifndef SPHEREGENERATION_H
#define SPHEREGENERATION_H

#include "HalfedgeDS.h"
#include "Eigen/Core"

/**
 * @author Luca Castelli Aleardi (2019)
 */
class SphereGeneration
{

public:
	/** 
	 * Initialize the data structures
	 **/
	SphereGeneration(Eigen::MatrixXd &V_original, Eigen::MatrixXi &F_original, HalfedgeDS &mesh);

    SphereGeneration() = default;
    
    ~SphereGeneration() = default;

	/** 
	 * Perform the subdivision of the mesh (just perform one round subdivision). <b>
	 * As result, the resulting subdivided mesh is stored in arrays 'V1' and 'F1' 
	 **/
	void subdivide();

	/** 
	 * Return the number of half-edges
	 **/
	Eigen::MatrixXd getVertexCoordinates()
	{
		return V1;
	}

	/** 
	 * Return the number of faces
	 **/
	Eigen::MatrixXi getFaces()
	{
		return F1;
	}

	/** 
	 * Print the combinatorial information of the subdivided mesh <b>
	 * verbosity=0: print only the number of vertices and faces <b>
	 * verbosity=1: print all incidence relations
	 **/
	void print(int verbosity);
	

private:

	/**
	 * Compute the midpoint of the given half-edge 'h=(u,v)'
	 */
	Eigen::MatrixXd computeEdgePoint(int h);

	/** 
     * Half-edge representation of the original input mesh 
     * */
	HalfedgeDS *he;

	Eigen::MatrixXd V; // vertex coordinates of the original input mesh

	/** faces/vertex incidence relations in the original mesh */
	Eigen::MatrixXi F; // REMARK: not needed if using the half-edge data structure

	int nVertices, nFaces; // number of vertices, faces in the new subdivided mesh
	Eigen::MatrixXd V1;		   // vertex coordinates of the new subdivided mesh
	Eigen::MatrixXi F1;		   // faces of the new subdivided mesh
};


#endif