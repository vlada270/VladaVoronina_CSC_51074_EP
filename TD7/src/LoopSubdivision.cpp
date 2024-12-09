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
    int f = F_original.rows();

    // TO BE COMPLETED (initialize arrays V1 and F1)
    nVertices = n + e; // Total vertices in V1
    nFaces = 4 * f;    // Total faces in F1
    V1.resize(nVertices, 3);
    F1.resize(nFaces, 3);
}


/** 
 * Perform the subdivision of the mesh (just perform one round subdivision). <b>
 * As result, the resulting subdivided mesh is stored in arrays 'V1' and 'F1' 
 **/

int getMidpointIndex(int v1, int v2, const std::map<std::pair<int, int>, int>& midPoints, int n) {
    if (v1 > v2) std::swap(v1, v2);
    auto it = midPoints.find(std::make_pair(v1, v2));
    return n + it->second;
}

void LoopSubdivision::subdivide()
{
    std::cout << "Performing one round of Loop subdivision" << std::endl;
    int e = he->sizeOfHalfedges() / 2; // Number of edges in the original mesh
    int n = he->sizeOfVertices();      // Number of vertices in the original mesh
    int F = he->sizeOfFaces();         // Number of faces in the original mesh

    //TODO: implement the subdivision algorithm
    V1.topRows(n) = V;

    std::map<std::pair<int, int>, int> mapEd;
    int count = 0;
    for (int i = 0; i < 2 * e; i++) {
        int v0 = he->getTarget(i);
        int v1 = he->getTarget(he->getOpposite(i));
        std::pair<int, int> edge = std::minmax(v0, v1);
        if (mapEd.find(edge) == mapEd.end())
        {
            V1.row(n + count) = computeEdgePoint(i);
            mapEd[edge] = count;
            count++;
        }
    }

    for (int i = 0; i < F; i++)
    {
        int h = he->getEdgeInFace(i);
        int h1 = he->getNext(h);
        int h2 = he->getNext(h1);

        int v1 = he->getTarget(h);
        int v2 = he->getTarget(h1);
        int v3 = he->getTarget(h2);

        int v4 = getMidpointIndex(v1, v2, mapEd, n);
        int v5 = getMidpointIndex(v2, v3, mapEd, n);
        int v6 = getMidpointIndex(v1, v3, mapEd, n);

        F1.row(4 * i) = Eigen::Vector3i(v1, v4, v6);
        F1.row(4 * i + 1) = Eigen::Vector3i(v4, v2, v5);
        F1.row(4 * i + 2) = Eigen::Vector3i(v4, v5, v6);
        F1.row(4 * i + 3) = Eigen::Vector3i(v6, v5, v3);
    }

    for (int i = 0; i < n; i++) {
        V1.row(i) = updateOriginalPoint(i);
    }
}

double LoopSubdivision::calculateGaussianCurvature(int v)
{
    double angleSum = 0.0;
    int edge = he->getEdge(v);
    int adjacentEdge = he->getNext(he->getOpposite(edge));

    do {
        int face = he->getFace(edge);
        double angle = 2 * M_PI / 3;
        angleSum += angle;

        edge = he->getNext(he->getOpposite(adjacentEdge));
        adjacentEdge = he->getNext(he->getOpposite(adjacentEdge));
    } while (adjacentEdge != edge);

    double curvature = 2 * M_PI - angleSum;
    return curvature;
}

void LoopSubdivision::subdivide_adaptive()
{
    std::cout << "Performing one round subdivision" << endl;
    int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
    int n = he->sizeOfVertices();      // number of vertices in the original mesh
    int F = he->sizeOfFaces();

    std::vector<double> curvature(n, 0.0);
    double threshold = 0.1;

    for (int v = 0; v < n; v++) {
        double K = calculateGaussianCurvature(v);
        curvature[v] = K;
    }

    V1.topRows(n) = V;

    std::map<std::pair<int, int>, int> mapEd;
    int count = 0;
    for (int i = 0; i < 2 * e; i++) {
        int v1 = he->getTarget(i);
        int v2 = he->getTarget(he->getOpposite(i));
        std::pair<int, int> edge = std::minmax(v1, v2);

        if (curvature[v1] > threshold || curvature[v2] > threshold) {
            if (mapEd.find(edge) == mapEd.end()) {
                V1.row(n + count) = computeEdgePoint(i);
                mapEd[edge] = count;
                count++;
            }
        }
    }

    for (int i = 0; i < F; i++) {
        int h = he->getEdgeInFace(i);
        int h1 = he->getNext(h);
        int h2 = he->getNext(h1);

        int v1 = he->getTarget(h);
        int v2 = he->getTarget(h1);
        int v3 = he->getTarget(h2);

        int v4 = getMidpointIndex(v1, v2, mapEd, n);
        int v5 = getMidpointIndex(v2, v3, mapEd, n);
        int v6 = getMidpointIndex(v1, v3, mapEd, n);

        F1.row(4 * i) = Eigen::Vector3i(v1, v4, v6);
        F1.row(4 * i + 1) = Eigen::Vector3i(v4, v2, v5);
        F1.row(4 * i + 2) = Eigen::Vector3i(v4, v5, v6);
        F1.row(4 * i + 3) = Eigen::Vector3i(v6, v5, v3);
    }

    for (int i = 0; i < n; i++) {
        V1.row(i) = updateOriginalPoint(i);
    }


    
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
    return 3.0 / 8 * V.row(he->getTarget(h)) + 3.0 / 8 * V.row(he->getTarget(he->getOpposite(h))) + 1.0 / 8 * V.row(he->getTarget(he->getNext(h))) + 1.0 / 8 * V.row(he->getTarget(he->getNext(he->getOpposite(h))));
}

/**
 * Given a vertex 'v' of the original mesh, compute and return its new coordinates
 */

double computeAlpha(int degree)
{
    if (degree == 3) {
        return 3.0 / 16;
    } else if (degree > 3) {
        return 3.0 / (8 * degree);
    }
    return 0.0;
}

MatrixXd LoopSubdivision::updateOriginalPoint(int v)
{
    //TODO
    MatrixXd newPoint = MatrixXd::Zero(1, 3);
    double alpha = computeAlpha(vertexDegree(v));
    int edge = he->getOpposite(he->getEdge(v));
    int adjacentEdge = he->getNext(he->getOpposite(edge));
    newPoint += alpha * V.row(he->getTarget(edge));
    newPoint += (1.0 - vertexDegree(v) * alpha) * V.row(v);
    while (adjacentEdge != edge)
    {
        int adjacentVertex = he->getTarget(adjacentEdge);
        newPoint += alpha * V.row(adjacentVertex);
        adjacentEdge = he->getNext(he->getOpposite(adjacentEdge));
    }

    return newPoint;
}

int LoopSubdivision::vertexDegree(int v)
{
    //TODO
    int count = 1;
    int edge = he->getOpposite(he->getEdge(v));
    int adjacentEdge = he->getNext(he->getOpposite(edge));
    while (adjacentEdge != edge) {
        adjacentEdge = he->getNext(he->getOpposite(adjacentEdge));
        count++;
    }
    return count;
}


