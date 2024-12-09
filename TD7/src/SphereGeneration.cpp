#include "SphereGeneration.h"
#include <iostream>
#include <map>

using Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXd;
using std::map, std::cout, std::endl;

SphereGeneration::SphereGeneration(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
{
    he = &mesh;
    V = V_original;
    F = F_original; // NOT NEEDED if using the half-edge data structure
    int e = he->sizeOfHalfedges() / 2; // number of edges in the original mesh
    int v = V_original.rows();         // number of vertices in the original mesh
    int f = F_original.rows();

    // Initialize arrays V1 and F1
    nVertices = v + e; // Total vertices in V1
    nFaces = 4 * f;    // Total faces in F1
    V1.resize(nVertices, 3);
    F1.resize(nFaces, 3);
}

void SphereGeneration::subdivide()
{
    cout << "Performing one round subdivision" << endl;
    int e = he->sizeOfHalfedges() / 2; // Number of edges in the original mesh
    int v = he->sizeOfVertices();      // Number of vertices in the original mesh
    int F = he->sizeOfFaces();         // Number of faces in the original mesh

    V1.topRows(v) = V;

    std::map<std::pair<int, int>, int> mapEd;  // Declare the map for edges to midpoint indices
    int i = v;
    int j = 0;

    for (int h = 0; h < he->sizeOfHalfedges(); h++)
    {
        int v0 = he->getTarget(h);
        int v1 = he->getTarget(he->getOpposite(h));
        std::pair<int, int> edge = std::minmax(v0, v1);
        if (mapEd.find(edge) == mapEd.end())
        {
            V1.row(i) = (V.row(v0) + V.row(v1)).normalized();
            mapEd[edge] = i;
            i++;
        }
    }
    for (int f = 0; f < F; f++)
    {
        int e0 = he->getEdgeInFace(f);
        int e1 = he->getNext(e0);
        int e2 = he->getNext(e1);
        int v0 = he->getTarget(e0);
        int v1 = he->getTarget(e1);
        int v2 = he->getTarget(e2);
        int m0 = mapEd[std::make_pair(std::min(v0, v1), std::max(v0, v1))];
        int m1 = mapEd[std::make_pair(std::min(v1, v2), std::max(v1, v2))];
        int m2 = mapEd[std::make_pair(std::min(v2, v0), std::max(v2, v0))];
        F1.row(j++) = Eigen::Vector3i(v0, m0, m2);
        F1.row(j++) = Eigen::Vector3i(m0, v1, m1);
        F1.row(j++) = Eigen::Vector3i(m1, v2, m2);
        F1.row(j++) = Eigen::Vector3i(m0, m1, m2);
    }
}

void SphereGeneration::print(int verbosity)
{
    cout << "\tn=" << nVertices << ", f=" << nFaces << endl;

    if (verbosity > 0) {
        // Print all vertex coordinates and face/vertex incidence relations
        for (int i = 0; i < nVertices; i++) {
            cout << "v" << i << ": " << V1.row(i) << endl;
        }

        cout << "new faces: " << nFaces << endl;
        for (int i = 0; i < nFaces; i++) {
            cout << "f" << i << ": " << F1.row(i) << endl;
        }
    }
}

Eigen::MatrixXd SphereGeneration::computeEdgePoint(int h)
{
    Eigen::RowVector3d midpoint = (V.row(he->getTarget(h)) + V.row(he->getTarget(he->getOpposite(h)))) / 2.0;
    midpoint.normalize();

    return midpoint;
}



