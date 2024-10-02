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
#ifndef LOOPSUBDIVISION_H
#define LOOPSUBDIVISION_H



#include "Eigen/Core"
#include "HalfedgeDS.h"
#include <memory>

using namespace Eigen;

/**
 * @author Luca Castelli Aleardi (2019)
 */
class LoopSubdivision
{

public:
    
    LoopSubdivision(Eigen::MatrixXd &V_original, Eigen::MatrixXi &F_original, HalfedgeDS &mesh);
    
    LoopSubdivision() = default;

    ~LoopSubdivision() = default;

    void subdivide();

    void subdivide_adaptive();
    
    MatrixXd getVertexCoordinates();

    MatrixXi getFaces();
    
    void print(int verbosity);

private:

    MatrixXd computeEdgePoint(int h);

    MatrixXd updateOriginalPoint(int v);

    int vertexDegree(int v);
    

    HalfedgeDS *he; // half-edge representation of the subdivided mesh
    MatrixXd V; // vertex coordinates of the original input mesh

	/** faces/vertex incidence relations in the original mesh */
	MatrixXi F; // REMARK: not needed if using the half-edge data structure

    int nVertices, nFaces; // number of vertices, faces in the new subdivided mesh
    MatrixXd V1;          // vertex coordinates of the new subdivided mesh
    MatrixXi F1;          // faces of the new subdivided mesh
};
#endif