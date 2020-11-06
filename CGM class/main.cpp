#include "distribution_function.h"

using namespace std;

int main() {
	const int X = 9;
	const int Y = 9;
	
	using namespace distribution_function_template_space;
	using dif = distribution_function_template_D2Q9<X, Y>;
	using dic = distribution_function_CGM_D2Q9<X, Y>;
	dif f;
	dic f_r,f_b;

	for (int i = 1; i <= X; i++) {
		for (int j = 1; j <= Y; j++) {
			printf("%8.4f", (dif::InM*dif::M)(i,j));
		}
		cout << endl;
	}
	f.single_phase_collison_SRT();
	vector<double, X> vec1(1);

	matrix<double, X, Y> m(f.S);
	printf("\n\n");
	for (int i = 1; i <= X; i++) {
		for (int j = 1; j <= Y; j++) {
			printf("%8.4f", m(i, j));
		}
		cout << endl;
	}
	
	f_r.single_phase_collison_MRT(f_r,f_b);
	cout<< 1/f_r.nu;
	return 0;
}