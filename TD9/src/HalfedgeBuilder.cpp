


#include <unordered_map>
#include <chrono>
#include <list>
#include <iostream>
#include "HalfedgeBuilder.h"

//using namespace Eigen;
using namespace std;

VertexPair::VertexPair(const int &v1, const int &v2)
{
  first = v1;
  second = v2;
}

bool VertexPair::operator==(const VertexPair &p) const
{
  if (first == p.first && second == p.second)
    return true;
  else if (first == p.second && second == p.first)
    return true;
  else
    return false;
}

class VertexPairHashFunction
{
public:
  // Use sum of lengths of first and last names
  // as hash function.
  size_t operator()(const VertexPair &p) const
  {
    return p.first * p.second;
  }
};



HalfedgeDS HalfedgeBuilder::createMesh(int nV, const MatrixXi &F)
{
  int nF = F.rows(); // number of triangle faces
  int nE;            // total number of half-edges: including both inner and boundary edges
  int nB;            // number of boundary edges
  int d = 3;         // the degree is 3 for all faces (in a triangle mesh)

  unordered_map<VertexPair, int, VertexPairHashFunction> insertedHalfedges; // allows to efficiently retrieve an inserted half-edge

  //cout << "Building an array-based implementation of an halfedge representation from a shared-vertex representation" << endl;
  auto start = std::chrono::high_resolution_clock::now(); // for measuring time performances

  // preliminary step: iterate over all faces to count inner and boundary edges (boundary edges are counted once)
  //cout << "Counting edges...";
  int innerEdges = 0;
  nE = 0;
  for (int i = 0; i < nF; i++)
  {
    int edge;
    for (int j = 0; j < d; j++)
    {
      MatrixXi f = F.row(i); // the incidence relations of the i-th face
      VertexPair pair0(f(j), f((j + 1) % d));
      //cout << "f" << i << ", " << counter << ": processing halfedge " << pair0;

      edge = i * d + j;
      if (insertedHalfedges.find(pair0) == insertedHalfedges.end())
      {
        insertedHalfedges[pair0] = edge;
        nE = nE + 2; // count all half-edges: each edge must be counted twice
        //cout << " inserted, " << insertedHalfedges.size() << endl;
      }
      else
        innerEdges++; // count inner edges: they are shared by an edge
    }
  }
  nB = (nE / 2) - innerEdges;
  //cout << "done" << endl;
  // cout << "\tn=" << nV << ", e=" << (nE / 2) << ", f=" << nF;
  // cout << ", \tboundary edges=" << nB << ", inner edges=" << innerEdges << endl;

  HalfedgeDS result(nV, nE); // the resulting triangle mesh

  int i;
  int edgeCounter = 0;

  //cout << "Setting the target vertex and the incident face...";
  for (i = 0; i < nF; i++)
  {
    // creating and storing the 'd' half-edges in a face
    for (int j = 0; j < d; j++)
    {
      // setting incident (target) vertex
      MatrixXi f = F.row(i); // the incidence relations of the i-th face
      result.setVertex(edgeCounter, f((j + 1) % d));

      // setting the face contaning the hal-edge
      result.setFace(edgeCounter, i);
      edgeCounter++;
    }
  }
  //cout << "done" << endl;

  // setting incidence relations between half-edges in the same face
  //cout << "Setting 'next' references...";
  int indNext; // index of next half-edge
  for (int j = 0; j < edgeCounter; j++)
  {
    if (j % 3 != 2) // only works for triangular meshes
      indNext = (j + 1);
    else
      indNext = (j - 2);

    result.setNext(j, indNext); // set the 'next' half-edge in the same face
    result.setPrev(indNext, j); // set the 'previous' half-edge in the same face
  }
  //cout << "done" << endl;

  // setting opposite half-edge
  insertedHalfedges.clear(); // reset the hash map (the map must be empty)
  int edge;
  int nInsertedHalfedges = 0;
  //cout << "Setting 'opposite' references...";
  int counter = 0;
  for (i = 0; i < nF; i++)
  {
    for (int j = 0; j < d; j++)
    {
      MatrixXi f = F.row(i); // the incidence relations of the i-th face
      VertexPair pair0(f(j), f((j + 1) % d));
      //cout << "f" << i << ", " << counter << ": processing halfedge " << pair0;

      edge = i * d + j;
      if (insertedHalfedges.find(pair0) == insertedHalfedges.end())
      {
        insertedHalfedges[pair0] = edge;
        //cout << " inserted, " << insertedHalfedges.size() << endl;
        nInsertedHalfedges++;
      }
      else
      {
        int eOpposite = insertedHalfedges[pair0];
        //cout << " not inserted, opposite " << eOpposite << endl;
        result.setOpposite(edge, eOpposite);
        result.setOpposite(eOpposite, edge);
      }

      counter++;
    }
  }
  //cout << "done (" << nInsertedHalfedges << " inserted edges)" << endl;

  // for vertices: set the incident half-edge (having the vertex as target)
  //cout << "Setting half-edges incident to edges...";
  int vertex, oppPrev;
  for (edge = 0; edge < nE; edge++)
  {
    vertex = result.getTarget(edge);
    if (vertex != -1)
      result.setEdge(vertex, edge);
  }
  //cout << "done" << endl;

  if (nB > 0) // process boundary edges if any
    addBoundaryEdges(result);

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  // std::cout << "Construction time: " << elapsed.count() << " s\n";

  return result;
}

/**
 * Add and set boundary half-edges to the representation
 * Remark: it allows to deal with meshes having boundaries
*/
void HalfedgeBuilder::addBoundaryEdges(HalfedgeDS mesh)
{
  int nE = mesh.sizeOfHalfedges(); // total number of halfedges the data structure must store
  std::list<int> boundaryEdges{};       // a list storing boundary edges to process
  int nInsertedHalfedges = 0;      // counts the number of inserted halfedges
  int bCount = 0;                  // count the number of boundaries
  int edgeCounter;                 // counter for enumerating half-edges

  for (int e = 0; e < nE; e++) // store all boundary edges
  {                            // Remark: boundary edges are those for which the 'opposite' is not defined, but the 'next' is already defines
    if (mesh.getOpposite(e) == -1 && mesh.getNext(e) != -1)
      boundaryEdges.push_back(e);
  }
  edgeCounter = nE - boundaryEdges.size();

  while (boundaryEdges.empty() == false) // set opposite references of all boundary edges
  {
    int firstEdge = boundaryEdges.front(); // get a boundary edge
    boundaryEdges.pop_front();             // remove the boundary edge

    int opFirstEdge = edgeCounter; // the index of the half-edge to be added
    edgeCounter++;                 // increment the counter (for adding the new half-edge)
    mesh.setOpposite(opFirstEdge, firstEdge);
    mesh.setOpposite(firstEdge, opFirstEdge);
    mesh.setVertex(opFirstEdge, mesh.getTarget(mesh.getNext((mesh.getNext(firstEdge)))));
    mesh.setFace(opFirstEdge, -1); // not required: just for the sake of clarity
  }

  for (int e = 0; e < nE; e++) // iterate over all edges
  {                            // Remark: boundary edges are those for which the 'opposite' is not defined, but the 'next' is already defines
    if (mesh.getFace(e) == -1) {
      int next=getNextBoundaryHalfedge(mesh, e);
      mesh.setNext(e, next);
      mesh.setPrev(next, e);
    }
  }
}

/**
 * Given a boundary exterior half-edge exterior 'e', returns the next boundary half-edge around the same boundary cycle
 * Warning: the mesh is supposed to be manifold (boundary cycles are disjoint)
 */
int HalfedgeBuilder::getNextBoundaryHalfedge(HalfedgeDS he, int e)
{
    if (he.getFace(e) != -1) // the starting edge must be a boundary (exterior) half-edge: its incidence face is not defined
      return -1;             // error

    // turn around the target vertex of 'e', starting from the half-edge pEdge
    int pEdge = he.getOpposite(e);
    while (he.getFace(pEdge) != -1)
    {
      //System.out.print(" pEdge:"+printEdge(pEdge));
      pEdge = he.getOpposite(he.getNext(he.getNext(pEdge)));
    }
    //System.out.println("next boundary edge:"+printEdge(pEdge));
    return pEdge;
}

