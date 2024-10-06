#include "Transform.h"
#include <iostream>


using Eigen::MatrixXd, Eigen::Vector3d, Eigen::RowVector3d, Eigen::Quaterniond;

void Transform::uniform_scale(double s){
	Eigen::MatrixXd m(4, 4);
    //TODO
    m <<    s, 0, 0, 0,
            0, s, 0, 0,
            0, 0, s, 0,
            0, 0, 0, 1;

    M = m;
}

void Transform::scale(double s1,double s2,double s3){
	Eigen::MatrixXd m(4, 4);
    //TODO
    m <<    s1, 0, 0, 0,
            0, s2, 0, 0,
            0, 0, s3, 0,
            0, 0, 0, 1;

    M = m;
}

void Transform::translate(Eigen::Vector3d t){
	Eigen::MatrixXd m(4, 4);
    //TODO
    m <<    1, 0, 0, t(0),
            0, 1, 0, t(1),
            0, 0, 1, t(2),
            0, 0, 0, 1;

    M = m;
}

void Transform::rotate_around_axis(double theta, int axis){
	Eigen::MatrixXd m(4, 4);
    //TODO
    m <<    1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

    int i1 = (axis + 1) % 3;
    int i2 = (axis + 2) % 3;

    m(i1, i1) = cos(theta);
    m(i2, i2) = cos(theta);
    m(i1, i2) = -sin(theta);
    m(i2, i1) = sin(theta);


    M = m;
}

Eigen::MatrixXd Transform::get_matrix(){
	return this->M;
}

void Transform::set_matrix(Eigen::MatrixXd m){
	this->M = m;
}


Transform Transform::operator*(const Transform& T){
	//TODO
    Transform transform;
    transform.M = this->M * T.M;
    return transform;
}

Transform Transform::operator=(const Transform& T){
	return(Transform(T.M));
}

Transform Transform::operator*(const Eigen::MatrixXd& M_in){
	//TODO
    Transform transform;
    transform.M = this->M * M_in;
    return transform;
}

Transform Transform::operator=(const Eigen::MatrixXd& M){
	//TODO
    this->M = M;
    return *this;
}

Transform Transform::compose(const Transform& T){
	//TODO
    return (*this) * T;
}

bool Transform::operator==(const Transform& T){
	//TODO
    return this->M.isApprox(T.M);  // Eigen has an isApprox function for approximate equality
}

void Transform::apply_transform(Eigen::MatrixXd& V){

	//TODO
    Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(V.rows(), 1);
    Eigen::MatrixXd V_homogeneous(V.rows(), 4);
    V_homogeneous << V, ones;

    V_homogeneous = (M * V_homogeneous.transpose()).transpose();
    V = V_homogeneous.leftCols(3);
}

void Transform::print_quaternion(Quaterniond q) {
	std::cout << "(" << q.w() << ", " << q.x() << ", "<< q.y() << ", "<< q.z() <<std::endl;
}

void Transform::rotate_with_quaternions(Eigen::MatrixXd& V, Eigen::Vector3d axis, double theta){
	// To be COMPLETED
    axis.normalize();
    Eigen::Quaterniond q(Eigen::AngleAxisd(theta, axis));
    for (int i = 0; i < V.rows(); ++i) {
        Eigen::Vector3d point = V.row(i);
        Eigen::Quaterniond v(0, point.x(), point.y(), point.z());
        Eigen::Quaterniond rotation = q * v * q.inverse();
        Eigen::Vector3d rotated_point(rotation.x(), rotation.y(), rotation.z());
        V.row(i) = rotated_point;
    }
}