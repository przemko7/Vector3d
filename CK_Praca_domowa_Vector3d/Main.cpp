#include <iostream>
#include "Vector3d.h"

int main() {
	Vector3d *A = new Vector3d;
	A->SetVector3d(2, 5, 10);
	A->Write();
	delete A;
	return 0;
}