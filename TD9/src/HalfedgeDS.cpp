

#include <iostream>
#include <ostream>
#include "HalfedgeDS.h"

using namespace std;


HalfedgeDS::HalfedgeDS(){};
/** 
 * Allocate the memory space for storing the data structure.
 * This version of the constructor does not allocate face indices (for reducing memory requirements)
 **/
HalfedgeDS::HalfedgeDS(int n, int h)
{
	nVertices = n;
	nHalfedges = h;
	T = new int[nHalfedges * sizeT];
	for (int i = 0; i < nHalfedges * sizeT; i++)
		T[i] = -1; // "-1" means that reference is NOT DEFINED (null)

	incidentEdge = new int[nVertices];
}

/** 
 * Set the opposite half-edge 
 **/
void HalfedgeDS::setOpposite(int e, int eOpposite)
{
	T[e * sizeT] = eOpposite;
}

/** 
 * Set the next half-edge: the half-edge 'eNext' following 'e' in the same face 
 **/
void HalfedgeDS::setNext(int e, int eNext)
{
	T[sizeT * e + 1] = eNext;
}

/** 
 * Set the (target) incident vertex of 'e' 
 **/
void HalfedgeDS::setVertex(int e, int v)
{
	T[sizeT * e + 2] = v;
}

/** 
 * Set the face containing the given half-edge
 **/
void HalfedgeDS::setFace(int e, int f)
{
	T[sizeT * e + 3] = f;
}

/** 
 * Set the face containing the given half-edge
 **/
void HalfedgeDS::setPrev(int e, int ePrev)
{
	T[sizeT * e + 4] = ePrev;
}

void HalfedgeDS::setEdge(int v, int e)
{
	incidentEdge[v] = e;
}

//--- methods for accessing data and navigating between mesh elements ---

/** 
 * Return the following edge in ccw orientation around the same face
 **/
int HalfedgeDS::getNext(int e)
{
	return T[sizeT * e + 1];
}

/** 
 * Return the opposite edge in the neighboring face (having opposite direction)
 **/
int HalfedgeDS::getOpposite(int e)
{
	return T[e * sizeT];
}

/** 
 * Return the previous edge in ccw orientation around the same face
 **/
int HalfedgeDS::getPrev(int e)
{
	return T[sizeT * e + 4];
}

/** 
 * Return the target vertex incident to the half-edge (its target vertex)
 **/
int HalfedgeDS::getTarget(int e)
{
	return T[sizeT * e + 2];
}

/** 
 * Return the face containing the half-edge
 **/
int HalfedgeDS::getFace(int e)
{
	return T[sizeT * e + 3];
}

/** 
 * Return a half edge incident to the vertex
 **/
int HalfedgeDS::getEdge(int v)
{
	return incidentEdge[v];
}

/** 
 * Return the number of vertices
 **/
int HalfedgeDS::sizeOfVertices()
{
	return nVertices;
}

/** 
 * Return the number of half-edges
 **/
int HalfedgeDS::sizeOfHalfedges()
{
	return nHalfedges;
}

/** 
 * Return the number of faces
 **/
int HalfedgeDS::sizeOfFaces()
{
	return nFaces;
}

/** 
 * Print the array T[] storing all references
 **/
void HalfedgeDS::print()
{
	for (int i = 0; i < nHalfedges; i++)
	{
		cout << "he" << i << ": \t" << T[sizeT*i] << "\t" << T[sizeT*i+1] << "\t" << T[sizeT*i+2] << "\t" << T[sizeT*i+3] << endl;
	}
}



