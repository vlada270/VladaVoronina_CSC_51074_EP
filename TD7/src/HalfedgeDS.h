#ifndef HALFEDGEDS_H
#define HALFEDGEDS_H

class HalfedgeDS
{

public:
	/** 
	 * Allocate the memory space for storing the data structure.
	 * This version of the constructor does not allocate face indices (for reducing memory requirements)
	 **/
	HalfedgeDS();
	
	HalfedgeDS(int n, int h);

	HalfedgeDS(int n, int h, int f);

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

	void setEdgeInFace(int f, int e);

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
	 * Return a half edge incident to the face
	 **/
	int getEdgeInFace(int f);

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

	int *faces; // for each face we store an incident halfedge
	const int sizeT=5;
};


#endif