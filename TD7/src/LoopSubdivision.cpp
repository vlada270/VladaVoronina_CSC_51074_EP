/* ========================================================================= *
 *                                                                           *
 *                       Luca Castelli Aleardi                       		 *
 *           Copyright (c) 2019, Ecole Polytechnique                		 *
 *           Department of Computer Science                  				 *
 *                          All rights reserved.                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of the course material developed for		             *
 *   INF574 Digital Representation and Analysis of Shapes (2019/20)			 *
 * ========================================================================= */
#include "igl/gaussian_curvature.h"
#include "LoopSubdivision.h"
#include <map>



using Eigen::MatrixXd, Eigen::MatrixXi;
using std::map, std::cout, std::endl;


/** 
 * Initialize the data structures
 **/
LoopSubdivision::LoopSubdivision(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
{
    he = &mesh;
    V = V_original;
    F = F_original; // NOT NEEDED if using the half-edge data structure
    int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
    int n = V_original.rows();         // number of vertices in the original mesh

    // TO BE COMPLETED (initialize arrays V1 and F1)
    

}


/** 
 * Perform the subdivision of the mesh (just perform one round subdivision). <b>
 * As result, the resulting subdivided mesh is stored in arrays 'V1' and 'F1' 
 **/
void LoopSubdivision::subdivide()
{
    std::cout << "Performing one round subdivision" << endl;
    int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
    int n = he->sizeOfVertices();      // number of vertices in the original mesh
    int F = he->sizeOfFaces();         // number of vertices in the original mesh
    //TODO: implement the subdivision algorithm

}

void LoopSubdivision::subdivide_adaptive()
{
    std::cout << "Performing one round subdivision" << endl;
    int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
    int n = he->sizeOfVertices();      // number of vertices in the original mesh

    //Perform adaptive subdivision


    
}	
    
/** 
 * Return the number of half-edges
 **/
MatrixXd LoopSubdivision::getVertexCoordinates()
{
    return V1;
}

/** 
 * Return the number of faces
 **/
MatrixXi LoopSubdivision::getFaces()
{
    return F1;
}

/** 
 * Print the combinatorial information of the subdivided mesh <b>
 * verbosity=0: print only the number of vertices and faces <b>
 * verbosity=1: print all incidence relations
 **/
void LoopSubdivision::print(int verbosity)
{
    cout << "\tn=" << nVertices << ", f=" << nFaces << endl;

    if (verbosity > 0) // print all vertex coordinates and face/vertex incidence relations
    {
        for (int i = 0; i < nVertices; i++)
        {
            cout << "v" << i << ": " << V1.row(i) << endl;
        }

        std::cout << "new faces: " << nFaces << endl;
        for (int i = 0; i < nFaces; i++)
        {
            cout << "f" << i << ": " << F1.row(i) << endl;
        }
    }
}

/**
 * Compute the midpoint of the given half-edge 'h=(u,v)'
 */
MatrixXd LoopSubdivision::computeEdgePoint(int h)
{
    //TODO
}

/**
 * Given a vertex 'v' of the original mesh, compute and return its new coordinates
 */
MatrixXd LoopSubdivision::updateOriginalPoint(int v)
{
    //TODO
}

int LoopSubdivision::vertexDegree(int v)
{
    //TODO
}


