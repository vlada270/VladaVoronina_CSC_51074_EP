#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "Eigen/Dense"

// In c++ classes can be defined in header files. This is a good practice if you have a small class that is only used in one file. If you have a class that is used in multiple files, you should define it in a header file and implement it in a cpp file. This is what we do here. We define the Transform class in Transform.h and implement it in Transform.cpp. This is a good practice because it makes your code more modular and easier to read. It also makes it easier to reuse your code in other projects.

//In c++ there are class and struct. The only difference between them is that in a class everything is private by default, while in a struct everything is public by default. In this case we use a struct, because we want everything to be public. This is a matter of personal preference, you can also use a class here.


struct Transform
{
private:
	/* data */
public:
	Eigen::MatrixXd M; //This is a 4x4 matrix that stores the transformation matrix. We will use this matrix to store the transformation matrix of the mesh. The matrix is initialized to the identity matrix in the constructor of the class.
	
	Transform(){
		Eigen::MatrixXd m(4, 4);
		m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
		m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
		m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
		m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

    	M = m;
	};

	Transform(Eigen::MatrixXd m){
		this->M = m;
	}

	~Transform() = default;


	void uniform_scale( double s);
	void scale( double s1,double s2, double s3);
	void translate(Eigen::Vector3d t);
	void rotate_around_axis(double angle, int axis);
	Eigen::MatrixXd get_matrix();
	void set_matrix(Eigen::MatrixXd m);
	
	Transform operator*(const Transform& T);
	Transform operator=(const Transform& T);
	Transform operator*(const Eigen::MatrixXd& M);
	Transform operator=(const Eigen::MatrixXd& M);
	bool operator==(const Transform& T);
	Transform compose(const Transform& T);


	void apply_transform(Eigen::MatrixXd& V);

	void rotate_with_quaternions(Eigen::MatrixXd& V, Eigen::Vector3d axis, double angle);

	void print_quaternion(Eigen::Quaterniond q);
};





#endif // TRANSFORM_H