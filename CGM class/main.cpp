#include "distribution_function.h"

int main() {
	
	using namespace distribution_function_template_space;

	distribution_function_template_D2Q9<172, 91> f(1, 0.1);

	f.setarea("area.txt");

	for (int step = 0; step <= 500; step++) {
		std::cout << step << std::endl;
		f.single_phase_collison_SRT();
		f.streaming();
		f.bounce_back_halfway();
		f.inlet_velocity_3(vector<double, 2>({ 0.01,0 }));
		f.oulet_outlflow_Neumann_2();
		f.density.calculate(f);
		f.velocity.calculate(f, f.density, f.area);
		f.velocity.detect(f.area);
		//f.detect();
	}

	f.out_file("test.vtk");

	return 0;
}