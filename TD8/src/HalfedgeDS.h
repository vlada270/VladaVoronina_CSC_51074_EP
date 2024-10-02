#ifndef HALFEDGEDS_H
#define HALFEDGEDS_H
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

/**
 * @author Luca Castelli Aleardi (2019)
 * Minimal array-based implementation of a the Half-edge data structure for representing polygonal meshes<br>
 * Features:<br>
 * -) References are encoded with integer values (32 bits), stores in an array of int <br>
 * -) vertex coordinates are not stored <br>
 * -) no implementation of dynamic updates <br>
 */

class HalfedgeDS
{

public:
	/** 
	 * Allocate the memory space for storing the data structure.
	 * This version of the constructor does not allocate face indices (for reducing memory requirements)
	 **/
	HalfedgeDS(int n, int h);

	/** 
	 * Set the opposite half-edge 
	 **/
	void setOpposite(int e, int eOpposite);

	/** 
	 * Set the next half-edge: the half-edge 'eNext' following 'e' in the same face 
	 **/
	void setNext(int e, int eNext);

	/** 
	 * Set the (target) incident vertex of 'e' 
	 **/
	void setVertex(int e, int v);

	/** 
	 * Set the face containing the given half-edge
	 **/
	void setFace(int e, int f);

	/** 
	 * Set the face containing the given half-edge
	 **/
	void setPrev(int e, int ePrev);

	void setEdge(int v, int e);

	//--- methods for accessing data and navigating between mesh elements ---

	/** 
	 * Return the following edge in ccw orientation around the same face
	 **/
	int getNext(int e);

	/** 
	 * Return the opposite edge in the neighboring face (having opposite direction)
	 **/
	int getOpposite(int e);

	/** 
	 * Return the previous edge in ccw orientation around the same face
	 **/
	int getPrev(int e);

	/** 
	 * Return the target vertex incident to the half-edge (its target vertex)
	 **/
	int getTarget(int e);

	/** 
	 * Return the face containing the half-edge
	 **/
	int getFace(int e);

	/** 
	 * Return a half edge incident to the vertex
	 **/
	int getEdge(int v);

	/** 
	 * Return the number of vertices
	 **/
	int sizeOfVertices();

	/** 
	 * Return the number of half-edges
	 **/
	int sizeOfHalfedges();

	/** 
	 * Return the number of faces
	 **/
	int sizeOfFaces();

	/** 
	 * Print the array T[] storing all references
	 **/
	void print();

private:
	int nVertices, nHalfedges, nFaces; // number of vertices, halfedges and faces in the mesh

	int *T;			   // a table for storing references between half-edges: each halfedge is represented with three integer references
	int *incidentEdge; // for each vertex we store an incident (ingoing) halfedge

	const int sizeT=5;
};


#endif