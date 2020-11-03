#include "distribution_function.h"
namespace distribution_function_base_space {

	//distribution_function_base

	distribution_function_base::distribution_function_base(int x_, int y_, int z_, int Q_)
		: x(x_), y(y_), z(z_), Q(Q_)
	{
		long xyzQ = x * y * z * Q;
		distribution_function_base_p = new double[xyzQ];
	}

	double& distribution_function_base::operator()(int i, int j, int k, int q) {
		// (i,j,k,q) <-- [i-1][j-1][k-1][q] <-- [q][i-1][j-1][k-1] <-- q*x*y*z + (i-1)*y*z + (j-1)*z + (k-1)
		int index = q * x * y * z + (i - 1) * y * z + (j - 1) * z + (k - 1);
		return *(distribution_function_base_p + index);
	}

	distribution_function_base::~distribution_function_base() {
		delete[]distribution_function_base_p;
	}

}

namespace distribution_function_space {

	//distribution_function_D2Q9

	distribution_function_D2Q9::distribution_function_D2Q9(int x_, int y_)
		:distribution_function_base(x_, y_, 1, 9),
		w({ 4.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0 })
	{
		c[0] = { 0,0 };
		c[1] = { 1,0 }; c[2] = { 0,1 }; c[3] = { -1,0 }; c[4] = { 0,-1 };
		c[5] = { 1,1 }; c[6] = { -1,1 }; c[7] = { -1,-1 }; c[8] = { 1,-1 };
	}

	double& distribution_function_D2Q9::operator()(int i, int j, int q) {
		return this->operator()(i, j, 1, q);
	}

	
	void distribution_function_D2Q9::shift(int q) {
		
		return;
	}

	void distribution_function_D2Q9::streaming() {

		for (int q = 1; q <= 8; q++) {
			switch (q)
			{
			case 1:
			case 2:
			case 3:
			case 4:
			case 5:
			case 6:
			case 7:
			case 8:
			default:
				break;
			}
		}

		return;
	}

}

namespace distribution_function_space {

	//velocity_field_2D

	velocity_field_2D::velocity_field_2D(int x_, int y_) :velocity_field<2>(x_,y_,1){}

	distribution_function_base_space::vector<double, 2>& velocity_field_2D::operator()(int i, int j) {
		return this->operator()(i, j, 1);
	}
}