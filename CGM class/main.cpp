#include "distribution_function.h"

using namespace std;

int main() {
	const int X = 7;
	const int Y = 5;

	using namespace distribution_function_template_space;
	using namespace distribution_function_template_space::vector2D_field_space;
	using namespace distribution_function_template_space::scalar_field_space;
	distribution_function_template_D2Q9<X, Y> f, f_r(1.0), f_b(0);
	density_field<X, Y> density(3), density_r(4), density_b(5);
	velocity2D_field<X, Y> velocity;
	force2D_field<X, Y> force(1.3, 2);
	vector<double, 3> vec, vec2({ 1, 2, 4 });
	phase_field<X, Y> phase_field1(1);

	

	area_field<areatype, X, Y> area(areatype::F);

	area(1, 2) = areatype::S;

	//for (int i = 1; i <= X; i++) {
	//	for (int j = 1; j <= Y; j++) {
	//		area(i, j) = (i - 1) * Y + j;
	//	}
	//}

	for (int i = 0; i <= X+1; i++) {
		for (int j = 0; j <= Y+1; j++) {
			printf("%3d", area(i, j));
		}
		cout << endl;
	}
	
	
	
	for (int i = 1; i <= X; i++) {
		for (int j = 1; j <= Y; j++) {
			cout << phase_field1.calculate(density_r, density_b)(i, j) << " ";
		}
		cout << endl;
	}
	cout << "\n\n";

	using v2 = vector<double, 2>;
	vector<double, 2> a[9]{v2({ 1, 2 }), v2({ 3, 4 })};

	cout << int(areatype::F) << endl;

	area_field<areatype, X, Y> area2(areatype::S);

	cout << int(area2(1, 1)) << endl;

	f.velocity.calculate(f,density);

	return 0;
}