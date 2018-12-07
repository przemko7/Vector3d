#include <iostream>
#include "Vector3d.h"

int main() {
	Vector3d *A = new Vector3d;
	A->SetVector3d(2, 5, 10);
	double a=A->Scalar(5, 6, 3);
	delete A;
	return 0;
}