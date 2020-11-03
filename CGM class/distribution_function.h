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

		using distribution_function_base_space::distribution_function_base::operator();
		
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

namespace distribution_function_template_base_space {

	//distribution_function_base_template

	template <int X, int Y, int Z, int Q>
	class distribution_function_base_template {
	public:

		distribution_function_base_template()
			:distribution_function_base_template_p(new double[X * Y * Z * Q]) {}

		~distribution_function_base_template() {
			delete[] distribution_function_base_template_p;
		}

		double& operator()(int i, int j, int k, int q) {
			// (i,j,k,q) <-- [i-1][j-1][k-1][q] <-- [q][i-1][j-1][k-1] <-- q * X * Y * Z + (i - 1) * Y * Z + (j - 1) * Z + (k - 1)
			int index = q * X * Y * Z + (i - 1) * Y * Z + (j - 1) * Z + (k - 1);
			return *(distribution_function_base_template_p + index);
		}

	protected:

		double* distribution_function_base_template_p;

	};

}

namespace distribution_function_template_space {

	//D2Q9 distribution_functions_template--------------------------

	template <int X,int Y>
	class distribution_function_template_D2Q9  {
	public:

		distribution_function_template_D2Q9()
			:w({ 4.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0 }),
			distribution_function_template_p (new double [X*Y*9])
		{
			c[0] = { 0,0 };
			c[1] = { 1,0 }; c[2] = { 0,1 }; c[3] = { -1,0 }; c[4] = { 0,-1 };
			c[5] = { 1,1 }; c[6] = { -1,1 }; c[7] = { -1,-1 }; c[8] = { 1,-1 };

			for (int r = 0; r < X * Y * 9; r++) {
				*(distribution_function_template_p + r) = 0;
			}
		}

		distribution_function_template_D2Q9(double rho_initial)
			:w({ 4.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0 }),
			distribution_function_template_p(new double[X * Y * 9])
		{
			c[0] = { 0,0 };
			c[1] = { 1,0 }; c[2] = { 0,1 }; c[3] = { -1,0 }; c[4] = { 0,-1 };
			c[5] = { 1,1 }; c[6] = { -1,1 }; c[7] = { -1,-1 }; c[8] = { 1,-1 };
			
			for (int q = 0; q <= 8; q++)
				for (int i = 1; i <= X; i++)
					for (int j = 1; j <= Y; j++)
						(*this)(i, j, q) = w(q) * rho_initial;
		}

		~distribution_function_template_D2Q9() {
			delete[]distribution_function_template_p;
		}

		double& operator()(int i, int j, int q) {
			// (i,j,q) <-- [i-1][j-1][q] <-- [q][i-1][j-1] <-- q * X * Y + (i - 1) * Y + (j - 1)
			int index = q * X * Y + (i - 1) * Y + (j - 1);
			return *(distribution_function_template_p + index);
		}

		//friend distribution_function_template_D2Q9<X, Y> operator+(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2);

		void blend(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2); //This function is used in the color gradient model to get the total distribution function.

	private:

		distribution_function_base_space::vector<double, 9> w;
		distribution_function_base_space::vector<double, 2> c[9];
		double* distribution_function_template_p;

	public:
		void streaming();

	};

	template <int X, int Y>
	void distribution_function_template_D2Q9<X, Y>::streaming() {

		for (int q = 1; q <= 8; q++) {
			switch (q)
			{
			case 1:
				for (int i = X; i >= 2; i--) for (int j = 1; j <= Y; j++)
					this->operator()(i, j, q) = this ->operator()(i - 1, j, q); break;
			case 2:
				for (int i = 1; i <= X; i++) for (int j = Y; j >= 2; j--)
					this->operator()(i, j, q) = this ->operator()(i, j - 1, q); break;
			case 3:
				for (int i = 1; i <= X - 1; i++) for (int j = 1; j <= Y; j++)
					this->operator()(i, j, q) = this ->operator()(i + 1, j, q); break;
			case 4:
				for (int i = 1; i <= X; i++) for (int j = 1; j <= Y - 1; j++)
					this->operator()(i, j, q) = this ->operator()(i, j + 1, q); break;
			case 5:
				for (int i = X; i >= 2; i--) for (int j = Y; j >= 2; j--)
					this->operator()(i, j, q) = this ->operator()(i - 1, j - 1, q); break;
			case 6:
				for (int i = 1; i <= X - 1; i++) for (int j = Y; j >= 2; j--)
					this->operator()(i, j, q) = this ->operator()(i + 1, j - 1, q); break;
			case 7:
				for (int i = 1; i <= X - 1; i++) for (int j = 1; j <= Y - 1; j++)
					this->operator()(i, j, q) = this ->operator()(i + 1, j + 1, q); break;
			case 8:
				for (int i = X; i >= 2; i--) for (int j = 1; j <= Y - 1; j++)
					this->operator()(i, j, q) = this ->operator()(i - 1, j + 1, q); break;
			default:
				break;
			}
		}

		return;
	}

	//template <int X,int Y>
	//distribution_function_template_D2Q9<X, Y> operator+(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2) {
	//}

	template <int X,int Y>
	void distribution_function_template_D2Q9<X, Y>::blend(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2) {
		for (int q = 0; q <= 8; q++)
			for (int i = 1; i <= X; i++)
				for (int j = 1; j <= Y; j++)
					(*this)(i, j, q) = f1(i, j, q) + f2(i, j, q);
		return;
	}

}

#endif // !_DISTRIBUTION_FUNCTION_H_
