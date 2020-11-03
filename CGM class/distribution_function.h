#ifndef _DISTRIBUTION_FUNCTION_H_
#define _DISTRIBUTION_FUNCTION_H_
namespace distribution_function_base_space {

	//distribution_function_base

	class distribution_function_base {
	public:
		distribution_function_base() = default;
		distribution_function_base(int x_, int y_, int z_, int Q_);
		~distribution_function_base();

		double& operator()(int i, int j, int k, int q);

	protected:
		int x;
		int y;
		int z;
		int Q;
		double* distribution_function_base_p;
	};
}

namespace distribution_function_base_space {

	//vector template--------------------------------------
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

		vector(const T(&a)[N]) : v(new T[N]) {
			for (int i = 0; i < N; i++)
				this->v[i] = a[i];
		}

		~vector() {
			delete[] v;
		}

		template <typename T, int N>
		friend T operator*(const vector<T, N>& a, const vector<T, N>& b);

		template <typename T, int N>
		friend vector<T, N> operator+(const vector<T, N>& a, const vector<T, N>& b);

		T& operator()(int i) {
			//start from 0
			return v[i];
		}

		vector<T, N>& operator=(const vector<T, N>& a) {
			for (int i = 0; i < N; i++)
				this->v[i] = a.v[i];
			return *this;
		}

		vector<T, N>& operator=(const T(&a)[N]) {
			for (int i = 0; i < N; i++)
				this->v[i] = a[i];
			return *this;
		}

		vector<T, N>& operator+=(vector<T, N>& a) {
			for (int i = 0; i < N; i++)
				this->v[i] += a.v[i];
			return *this;
		}

	private:
		T* v;
	};

	// * for vecotr template
	template <typename T, int N>
	T operator*(const vector<T, N>& a, const vector<T, N>& b) {
		T temp(0);
		for (int i = 0; i < N; i++)
			temp += (a.v[i] * b.v[i]);
		return temp;
	}

	// + for vector template 
	template <typename T, int N>
	vector<T, N> operator+(const vector<T, N>& a, const vector<T, N>& b) {
		vector<T, N> temp_vector;
		for (int i = 0; i < N; i++)
			temp_vector.v[i] = a.v[i] + b.v[i];
		return temp_vector;
	}

	//velocity_field template------------------------------
	template <int N>
	class velocity_field {
	public:
		velocity_field() = default;
		velocity_field(int x_, int y_, int z_);
		~velocity_field();

		vector<double, N>& operator()(int i, int j, int k);

	private:
		int x;
		int y;
		int z;
		vector<double, N>* velocity_p;
	};

	template <int N>
	velocity_field<N>::velocity_field(int x_, int y_, int z_)
		:x(x_), y(y_), z(z_)
	{
		long xyz = x * y * z;
		velocity_p = new vector<double, N>[xyz];
	}

	template <int N>
	vector<double, N>& velocity_field<N>::operator()(int i, int j, int k) {
		// (i,j,k) <-- [i-1][j-1][k-1] <-- (i-1)*y*z + (j-1)*z + (k-1)
		int index = (i - 1) * y * z + (j - 1) * z + (k - 1);
		return *(velocity_p + index);
	}

	template <int N>
	velocity_field<N>::~velocity_field() {
		delete[]velocity_p;
	}

}

namespace distribution_function_space {

	//D2Q9 distribution_functions--------------------------

	class distribution_function_D2Q9 :public distribution_function_base_space::distribution_function_base {
	public:

		using distribution_function_base::operator();
		
		distribution_function_D2Q9() = default;
		distribution_function_D2Q9(int x_, int y_);

		double& operator()(int i, int j,int q);
		//friend distribution_function_D2Q9 operator+(distribution_function_D2Q9&f1, distribution_function_D2Q9&f2);
	private:

		distribution_function_base_space::vector<double, 9> w;
		distribution_function_base_space::vector<double, 2> c[9];

	public:
		void streaming();

	};

	//velocity_field_2D------------------------------------

	class velocity_field_2D: distribution_function_base_space::velocity_field<2>{
	public:

		using velocity_field::operator();

		velocity_field_2D() = default;
		velocity_field_2D(int x_, int y_);

		distribution_function_base_space::vector<double, 2>& operator()(int i, int j);

	};

}

#endif // !_DISTRIBUTION_FUNCTION_H_
