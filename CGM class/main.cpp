#include "distribution_function.h"
#define X 5
#define Y 5
using namespace std;

int main() {

	using namespace distribution_function_template_space;
	using namespace distribution_function_template_space::velocity2D_template_space;
	using namespace distribution_function_template_space::scalar_field_space;
	
	velocity2D_template<X, Y> velocity(0.1,133);
	scalar_field<X, Y> density(1);
	distribution_function_template_D2Q9<X, Y> f(1),feq(1);
	
	cout << velocity(1, 2)(0) << " " << velocity(1, 2)(1) << endl;
	cout << f(1, 2, 0) << endl;

	cout << "\n\n";
	for (int q = 0; q <= 8; q++) {
		printf("%7.5f ", feq(1, 2, q));
	}

	cout << "\n\n";
	for (int q = 0; q <= 8; q++) {
		printf("%7.5f ", feq.equilibrium(velocity, density)(1, 2, q));
	}

	feq.detect();

	return 0;
}