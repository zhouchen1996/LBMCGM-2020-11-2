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

	//define a vector for convenience
	template <typename T, int N>
	class vector {
	public:
		vector() = default;

		template <typename T, int N>
		friend T operator*(vector<T, N>& a, vector<T, N>& b);

		T& operator()(int i) {
			//start from 1
			return v[i - 1];
		}
	private:
		T v[N];
	};

	//overload multiplication related to template class vector
	template <typename T, int N>
	T operator*(vector<T, N>& a, vector<T, N>& b) {
		T temp(0);
		for (int i = 1; i <= N; i++)
			temp += (a(i) * b(i));
		return temp;
	}

	//D2Q9 distribution_function
	class distribution_function_D2Q9 :public distribution_function_base_space::distribution_function_base {
	public:

		using distribution_function_base::operator();

		template <typename T, int N>
		friend class vector;

		distribution_function_D2Q9() = default;
		distribution_function_D2Q9(int x_, int y_);

		double& operator()(int i, int j,int q);
		//void streaming();
	private:
		double w[9];
		vector<double, 2> c[9];
		vector<double, 2> velocity;
	};

	
}

#endif // !_DISTRIBUTION_FUNCTION_H_
