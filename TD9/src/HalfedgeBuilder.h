#ifndef HALFEDGEBUILDER_H
#define HALFEDGEBUILDER_H

#include "HalfedgeDS.h"
#include "Eigen/Dense"

using namespace Eigen;
//using namespace std;

struct VertexPair
{
  int first, second;

  VertexPair(const int &v1, const int &v2);

  bool operator==(const VertexPair &p) const;
  

};





class HalfedgeBuilder
{

  public:
      HalfedgeBuilder(){};
      ~HalfedgeBuilder() = default;

      /**
      * Efficient construction of the Half-edge DS for triangle meshes <br>
      * <br>
      * Remarks: <br>
      * -) it works only for triangle meshes WITHOUT boundaries <br>
      * -) it does not store the geometry (point locations) but only the combinatorial structure
      */
      HalfedgeDS createMesh(int nV, const MatrixXi &F);

      int getNextBoundaryHalfedge(HalfedgeDS he, int e);


  private:
      /**
       * Add and set boundary half-edges to the representation
       * Remark: it allows to deal with meshes having boundaries
      */
      void addBoundaryEdges(HalfedgeDS mesh);
  
};

#endif