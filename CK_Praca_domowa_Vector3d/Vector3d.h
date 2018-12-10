#pragma once
#include<vector>
using namespace std;

class Vector3d {
	public:
		Vector3d();
		Vector3d(double x, double y, double z);
		~Vector3d();

		void ErrorVectorLenghtCorrection(vector<double>& vec);

		double GetVector3dX();
		double GetVector3dY();
		double GetVector3dZ();
		vector<double> GetVector3dVector();

		void SetVector3d(double x, double y, double z);
		void SetVector3d(Vector3d &B);
		void SetVector3d(std::vector<double> &vector);

		void AddVector3d(double x, double y, double z);
		void AddVector3d(Vector3d &B);
		void AddVector3d(vector<double> &vector);

		Vector3d AddVectorsVector3d(vector<double> &vec1, vector<double> &vec2);
		Vector3d AddVectorsVector3d(vector<double> &vec1, double x, double y, double z);
		Vector3d AddVectorsVector3d(double x1, double y1, double z1, double x2, double y2, double z2);
		vector<double> AddVectorsVector(vector<double> &vec1, vector<double> &vec2);
		vector<double> AddVectorsVector(vector<double> &vec1, double x, double y, double z);
		vector<double> AddVectorsVector(double x1, double y1, double z1, double x2, double y2, double z2);

		double Scalar(double x, double y, double z);
		double Scalar(Vector3d &A);
		double Scalar(vector<double> &vec1);
		double Scalar(vector<double> &vec1, vector<double> &vec2);
		double Scalar(double x1, double y1, double z1, double x2, double y2, double z2);

		void VectorProduct(double x, double y, double z);
		void VectorProduct(Vector3d &B);
		void VectorProduct(std::vector<double> &vec1);

		void VectorProduct(double x1, double y1, double z1, double x2, double y2, double z2);
		void VectorProduct(Vector3d &A, Vector3d &B);
		void VectorProduct(vector<double> &vec1, vector<double> &vec2);

		void Write();

	private:
		double x_;
		double y_;
		double z_;

};