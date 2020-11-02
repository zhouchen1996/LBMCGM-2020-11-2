#ifndef _DISTRIBUTION_FUNCTION_H_
#define _DISTRIBUTION_FUNCTION_H_
namespace distribution_function_base_space {

	class distribution_function_base {
	public:
		distribution_function_base() = default;
		distribution_function_base(int x_, int y_, int z_, int Q_);
		~distribution_function_base();

		double& operator()(int i, int j, int k, int q);

	private:
		int x;
		int y;
		int z;
		int Q;
		double* distribution_function_base_p;
	};

}

namespace distribution_function_space {

	class distribution_function_D2Q9 :public distribution_function_base_space::distribution_function_base {
	public:
		distribution_function_D2Q9() = default;
		distribution_function_D2Q9(int x_, int y_);


	private:
		int x;
		int y;
	};
}

#endif // !_DISTRIBUTION_FUNCTION_H_
