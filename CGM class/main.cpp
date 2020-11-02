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

	distribution_function_space::vector<double, 3> ux, uy;
	ux(1) = 6; uy(1) = 3;
	ux(2) = 5; uy(2) = 1;
	ux(3) = 9; uy(3) = 4;
	
	cout<< ux * uy << endl;
	cout << 6 * 3 + 5 * 1 + 9 * 4 << endl;

	distribution_function_space::vector<double, 3> u3(ux);
	cout << u3(1) << " " << u3(2) << " " << u3(3) << endl;

	u3 = uy + ux;
	cout << ux(1) << " " << ux(2) << " " << ux(3) << endl;
	cout << uy(1) << " " << uy(2) << " " << uy(3) << endl;
	cout << u3(1) << " " << u3(2) << " " << u3(3) << endl;

	return 0;
}