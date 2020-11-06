#include "distribution_function.h"

using namespace std;

int main() {
	const int X = 5;
	const int Y = 5;

	using namespace distribution_function_template_space;
	using namespace distribution_function_template_space::vector2D_field_space;
	using namespace distribution_function_template_space::scalar_field_space;
	distribution_function_template_D2Q9<X, Y> f, f_r(1.0), f_b(0);
	density_field<X, Y> density(3),density_r(4),density_b(5);
	velocity2D_field<X, Y> velocity;
	force2D_field<X, Y> force(1.3,2);
	vector<double, 3> vec, vec2({1, 2, 4});
	phase_field<X, Y> phase_field1(1);
	area_field<char, X, Y> area('a');
	for (int i = 1; i <= X; i++) {
		for (int j = 1; j <= Y; j++) {
			cout <<phase_field1.calculate(density_r,density_b)(i, j) << " ";
		}
		cout << endl;
	}
	cout << "\n\n";

	return 0;
}