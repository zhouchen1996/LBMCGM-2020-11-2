#include "distribution_function.h"
namespace distribution_function_base_space {

	//base class distribution_function_base

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

	//D2Q9
	distribution_function_D2Q9::distribution_function_D2Q9(int x_, int y_)
		:distribution_function_base(x_, y_, 1, 9), x(x_), y(y_) {}

}