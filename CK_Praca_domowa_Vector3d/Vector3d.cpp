#include "Vector3d.h"
#include <iostream>

using namespace std;

Vector3d::Vector3d() {
	this->x_ = 0;
	this->y_ = 0;
	this->z_ = 0;
}

Vector3d::Vector3d(double x, double y, double z) {
	this->x_ = x;
	this->y_ = y;
	this->z_ = z;
}

Vector3d::~Vector3d() {

}

void Vector3d::ErrorVectorLenghtCorrection(vector<double>& vec)
{
	cerr << "Podany wektor nie mia? 3 elementów" << endl;
	cerr << "Zostana usuniete elemnety o indeksie wiekszym niz 2 lub dopisane zostana 0" << endl;

	while (vec.size() < 3) {
		vec.push_back(0);
	}

	if (vec.size() > 3) {
		double t_a, t_b, t_c;

		t_a = vec[0];
		t_b = vec[1];
		t_c = vec[2];

		vec.clear();

		vec.push_back(t_a);
		vec.push_back(t_b);
		vec.push_back(t_c);
	}
}

double Vector3d::GetVector3dX() {
	return this->x_;
}

double Vector3d::GetVector3dY() {
	return this->y_;
}

double Vector3d::GetVector3dZ() {
	return this->z_;
}

vector<double> Vector3d::GetVector3dVector()
{
	vector<double> vec(3);
	vec[0] = this->x_;
	vec[1] = this->y_;
	vec[2] = this->z_;
	return vec;
}

void Vector3d::SetVector3d(double x, double y, double z) {
	this->x_ = x;
	this->y_ = y;
	this->z_ = z;
}

void Vector3d::SetVector3d(Vector3d &B) {
	this->x_ = B.GetVector3dX();
	this->y_ = B.GetVector3dY();
	this->z_ = B.GetVector3dZ();
}

void Vector3d::SetVector3d(std::vector<double>& vector) {
	try {
		if (vector.size() != 3) {
			throw vector;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}
	this->x_ = vector[0];
	this->y_ = vector[1];
	this->z_ = vector[2];
}

void Vector3d::AddVector3d(double x, double y, double z) {
	this->x_ += x;
	this->y_ += y;
	this->z_ += z;
}

void Vector3d::AddVector3d(Vector3d &B) {
	this->x_ += B.GetVector3dX();
	this->y_ += B.GetVector3dY();
	this->z_ += B.GetVector3dZ();
}

void Vector3d::AddVector3d(vector<double>& vector)
{
	try {
		if (vector.size() != 3) {
			throw vector;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}
	this->x_ += vector[0];
	this->y_ += vector[1];
	this->z_ += vector[2];
}

Vector3d Vector3d::AddVectorsVector3d(vector<double>& vec1, vector<double>& vec2)
{
	try {
		if (vec1.size() != 3) {
			throw vec1;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}
	
	try {
		if (vec2.size() != 3) {
			throw vec2;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}

	Vector3d A;
	double x, y, z;
	x = vec1[0] + vec2[0];
	y = vec1[1] + vec2[1];
	z = vec1[2] + vec2[2];
	A.SetVector3d(x, y, z);
	return A;
}

Vector3d Vector3d::AddVectorsVector3d(vector<double>& vec1, double x, double y, double z)
{
	try {
		if (vec1.size() != 3) {
			throw vec1;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}

	Vector3d A;
	x += vec1[0];
	y += vec1[1];
	z += vec1[2];
	A.SetVector3d(x, y, z);
	return A;
}

Vector3d Vector3d::AddVectorsVector3d(double x1, double y1, double z1, double x2, double y2, double z2)
{
	Vector3d A;
	x1 += x2;
	y1 += y2;
	z1 += z2;
	A.SetVector3d(x1, y1, z1);
	return A;
}

vector<double> Vector3d::AddVectorsVector(std::vector<double>& vec1, std::vector<double>& vec2)
{
	try {
		if (vec1.size() != 3) {
			throw vec1;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}
	try {
		if (vec2.size() != 3) {
			throw vec2;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}

	vec1[0] += vec2[0];
	vec1[1] += vec2[1];
	vec1[2] += vec2[2];
	return vec1;
}

vector<double> Vector3d::AddVectorsVector(vector<double>& vec1, double x, double y, double z)
{
	try {
		if (vec1.size() != 3) {
			throw vec1;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}

	vec1[0] += x;
	vec1[1] += y;
	vec1[2] += z;
	return vec1;
}

vector<double> Vector3d::AddVectorsVector(double x1, double y1, double z1, double x2, double y2, double z2)
{
	vector<double> vec(3);
	vec[0] = x1 + x2;
	vec[1] = y1 + y2;
	vec[2] = z1 + z2;
	return vec;
}

double Vector3d::Scalar(double x, double y, double z) {
	double a, b, c;
	a = this->x_*x;
	b = this->y_*y;
	c = this->z_*z;
	return a + b + c;
}

double Vector3d::Scalar(Vector3d &A) {
	double a, b, c;
	a = A.GetVector3dX();
	b = A.GetVector3dY();
	c = A.GetVector3dZ();
	a *= this->x_;
	b *= this->y_;
	c *= this->z_;
	return a + b + c;
}

double Vector3d::Scalar(vector<double>& vec1) {
	try {
		if (vec1.size() != 3) {
			throw vec1;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}

	double a, b, c;
	a = this->x_*vec1[0];
	b = this->y_*vec1[1];
	c = this->z_*vec1[2];
	return a + b + c;
}

double Vector3d::Scalar(vector<double>& vec1, vector<double>& vec2) {
	try {
		if (vec1.size() != 3) {
			throw vec1;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}
	
	try {
		if (vec2.size() != 3) {
			throw vec2;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}

	double a, b, c;
	a = vec1[0] *vec2[0];
	b = vec1[1] *vec2[1];
	c = vec1[2] *vec2[2];
	return a + b + c;
}

double Vector3d::Scalar(double x1, double y1, double z1, double x2, double y2, double z2) {
	double a, b, c;
	a = x1 * x2;
	b = y1 * y2;
	c = z1 * z2;
	return a + b + c;
}

void Vector3d::VectorProduct(double x, double y, double z) {
	double a, b, c;
	a = this->y_*z - this->z_*y;
	b = -(this->x_*z - this->z_*x);
	c = this->x_*y - this->y_*x;
	this->SetVector3d(a, b, c);
}

void Vector3d::VectorProduct(Vector3d &B) {
	double a, b, c;

	a = this->y_*B.GetVector3dZ() - 
		this->z_*B.GetVector3dY();

	b = -(this->x_*B.GetVector3dZ() -
		this->z_*B.GetVector3dX());

	c = this->x_*B.GetVector3dY() -
		this->y_*B.GetVector3dX();

	this->SetVector3d(a, b, c);
}

void Vector3d::VectorProduct(std::vector<double>& vec1) {
	try {
		if (vec1.size() != 3) {
			throw vec1;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}

	double a, b, c;

	a = this->y_*vec1[2] -
		this->z_*vec1[1];

	b = -(this->x_*vec1[2] -
		this->z_*vec1[0]);

	c = this->x_*vec1[1] -
		this->y_*vec1[0];

	this->SetVector3d(a, b, c);
}

void Vector3d::VectorProduct(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double a, b, c;

	a = y1*z2 - 
		z1*y2;
	b = -(x1*z2 - 
		z1*x2);
	c = x1*y2 - 
		y1*x2;

	this->SetVector3d(a, b, c);
}

void Vector3d::VectorProduct(Vector3d & A, Vector3d & B)
{
	double a, b, c;

	a = A.GetVector3dY() * B.GetVector3dZ() -
		A.GetVector3dZ() * B.GetVector3dY();

	b = -(A.GetVector3dX() * B.GetVector3dZ() -
		A.GetVector3dZ() * B.GetVector3dX());

	c = A.GetVector3dX() * B.GetVector3dY() -
		A.GetVector3dY() * B.GetVector3dX();

	this->SetVector3d(a, b, c);
}

void Vector3d::VectorProduct(vector<double>& vec1, vector<double>& vec2)
{
	try {
		if (vec1.size() != 3) {
			throw vec1;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}

	try {
		if (vec2.size() != 3) {
			throw vec2;
		}
	}
	catch (std::vector<double> vec) {
		this->ErrorVectorLenghtCorrection(vec);
	}

	double a, b, c;
	a = vec1[1] *vec2[2] -
		vec1[2] *vec2[1];

	b = -(vec1[0] *vec2[2] -
		vec1[2] *vec2[0]);

	c = vec1[0] *vec2[1] -
		vec1[1] *vec2[0];

	this->SetVector3d(a, b, c);
}

void Vector3d::Write()
{
	printf("%.2f %.2f %.2f \n", this->x_, this->y_, this->z_);
}

