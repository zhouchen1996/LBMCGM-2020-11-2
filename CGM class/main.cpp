#include <iostream>
#include "distribution_function.h"

using namespace std;
int main() {
	distribution_function_space::distribution_function_D2Q9 f(3, 3);
	int temp = 0;
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int q = 0; q < 9; q++) {
				temp++;
				f(i, j, q) = temp;
			}
		}
	}
	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			for (int q = 0; q < 9; q++) {
				cout << f(i, j, q) << endl;
			}
		}
	}

	distribution_function_base_space::vector<double, 3> ux, uy;
	ux(0) = 6; uy(0) = 3;
	ux(1) = 5; uy(1) = 1;
	ux(2) = 9; uy(2) = 4;
	
	cout<< ux * uy << endl;
	cout << 6 * 3 + 5 * 1 + 9 * 4 << endl;

	distribution_function_base_space::vector<double, 3> u3(ux);
	cout << u3(0) << " " << u3(1) << " " << u3(2) << endl;

	u3 = uy + ux;
	cout << ux(0) << " " << ux(1) << " " << ux(2) << endl;
	cout << uy(0) << " " << uy(1) << " " << uy(2) << endl;
	cout << u3(0) << " " << u3(1) << " " << u3(2) << endl;
	cout << ux * uy + u3 * uy + 3 << endl;
	double aaa = (ux + uy) * u3;
	cout << aaa << endl;

	distribution_function_base_space::vector<double, 9> u4({55,32,3,32,45,32,543});
	cout << u4(0) << " " << u4(4) << " " << u4(8) << endl;

	u4 = { 0,1,2,3,4,5,6,7,8 };
	cout << u4(0) << " " << u4(4) << " " << u4(8) << endl;

	distribution_function_base_space::velocity_field<3> v(3, 3, 3);
	v(1, 1, 1) = { 1,2,1 };
	v(2, 3, 1) = { 1,1,1 };
	cout << v(1, 1, 1) * v(2,3,1) << endl;

	v(2, 3, 1) += v(1, 1, 1);
	cout << v(2, 3, 1)(0) << " " << v(2, 3, 1)(1) << " " << v(2, 3, 1)(2) << endl;
	v(2, 2, 2) = { 1,0,1 };
	cout << v(2, 2, 2) * v(2, 3, 1) << endl;

	distribution_function_space::velocity_field_2D v2d(4, 4);
	v2d(1, 1) = { 1,1 };
	v2d(2, 2) = { 2,2 };
	cout << v2d(1, 1) * v2d(2, 2) << endl;
	v2d(1, 1) += v2d(2, 2);
	cout << v2d(1, 1)(0) << " " << v2d(1, 1)(1) << endl;
	


	return 0;
}