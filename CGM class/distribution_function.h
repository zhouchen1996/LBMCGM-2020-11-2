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

		vector() 
			:v(new T[N])
		{
			for (int i = 0; i < N; i++) {
				v[i] = 0;
			}
		}

		vector(const vector<T, N>& a) : v(new T[N]) {
			for (int i = 0; i < N; i++)
				this->v[i] = a.v[i];
		}

		~vector() {
			delete[] v;
		}

		template <typename T, int N>
		friend T operator*(const vector<T, N>& a,const vector<T, N>& b);

		template <typename T, int N>
		friend vector<T, N> operator+(const vector<T, N>& a, const vector<T, N>& b);

		T& operator()(int i) {
			//start from 1
			return v[i - 1];
		}

		vector<T,N>& operator=(const vector<T, N>& a) {
			for (int i = 0; i < N; i++)
				this->v[i] = a.v[i];
			return *this;
		}

		vector<T, N>& operator+=(vector<T, N>& a) {
			for (int i = 0; i < N; i++)
				this->v[i] += a.v[i];
			return *this;
		}

	private:
		T *v;
	};

	//overload multiplication related to template class vector
	template <typename T, int N>
	T operator*(const vector<T, N>& a, const vector<T, N>& b) {
		T temp(0);
		for (int i = 0; i < N; i++)
			temp += (a.v[i] * b.v[i]);
		return temp;
	}

	template <typename T, int N>
	vector<T, N> operator+(const vector<T, N>& a,const vector<T, N>& b) {
		vector<T, N> temp_vector;
		for (int i = 0; i < N; i++)
			temp_vector.v[i] = a.v[i] + b.v[i];
		return temp_vector;
	}

	//D2Q9 distribution_functions
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
		vector<double,9> w;
		vector<double, 2> c[9];
	};

	class velocity {
	public:
		
	};

}

#endif // !_DISTRIBUTION_FUNCTION_H_
