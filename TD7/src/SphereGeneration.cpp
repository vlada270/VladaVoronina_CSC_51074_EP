#include "SphereGeneration.h"
#include <iostream>
#include <map>

using Eigen::MatrixXd,Eigen::MatrixXi, Eigen::VectorXd;
using std::map, std::cout, std::endl;

SphereGeneration::SphereGeneration(MatrixXd &V_original, MatrixXi &F_original, HalfedgeDS &mesh)
{
	//TODO
}

void SphereGeneration::subdivide()
{
	std::cout << "Performing one round subdivision" << endl;
	//TODO


}

void SphereGeneration::print(int verbosity)
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

MatrixXd SphereGeneration::computeEdgePoint(int h)
{
	//TODO 
}


