#include "distribution_function.h"
#define X 5
#define Y 5
using namespace std;

int main() {

	using namespace distribution_function_template_space;
	using namespace distribution_function_template_space::velocity2D_template_space;
	using namespace distribution_function_template_space::scalar_field_space;
	
	velocity2D_template<X, Y> velocity(0.1,0.4);
	distribution_function_template_D2Q9<X, Y> feq(1);
	density_field<X, Y> density(0.1);


	cout << "\n\n";
	for (int q = 0; q <= 8; q++) {
		printf("%7.5f ", feq(1, 2, q));
	}

	cout << "\n\n";
	for (int q = 0; q <= 8; q++) {
		printf("%7.5f ", feq.equilibrium(velocity, density)(1, 5, q));
	}

	feq.detect();

	density.calculate(feq);

	vector<double, 2> v({ 1,2 });
	cout << v(3) << endl;

	velocity2D_template<8, 8> veloc({ 1,0 });
	cout<<veloc(2, 0)(2);

	density(2, 6);

	return 0;
}