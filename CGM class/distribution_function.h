#ifndef _DISTRIBUTION_FUNCTION_H_
#define _DISTRIBUTION_FUNCTION_H_

//vector template
namespace distribution_function_template_space {
	//important!!!!!!!!
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

}

//velocity2D_template
namespace distribution_function_template_space {

	namespace velocity2D_template_space {

		template <int X, int Y>
		class velocity2D_template {
		public:

			velocity2D_template(double initial_velocity_x = 0, double initial_velocity_y = 0) :velocity2D_template_p(new distribution_function_template_space::vector<double, 2>[X * Y]) {
				for (int r = 0; r < X * Y; r++)
					*(velocity2D_template_p + r) = { initial_velocity_x,initial_velocity_y };
			}

			~velocity2D_template() {
				delete[]velocity2D_template_p;
			}

			distribution_function_template_space::vector<double, 2>& operator()(int i, int j) {
				// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
				int index = (i - 1) * Y + (j - 1);
				return *(velocity2D_template_p + index);
			}

		protected:

			distribution_function_template_space::vector<double, 2>* velocity2D_template_p;

		};

	}

}

//scalar_field  template
namespace distribution_function_template_space {

	namespace scalar_field_space {

		template<int X,int Y>
		class scalar_field {
		public:

			scalar_field(double scalar_initial = 0)
				: scalar_p(new double[X*Y])
			{
				for (int r = 0; r < X * Y; r++)
					*(scalar_p + r) = scalar_initial;
			}

			~scalar_field(){
				delete[]scalar_p;
			}

			double& operator()(int i,int j) {
				// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
				int index = (i - 1) * Y + (j - 1);
				return *(scalar_p + index);
			}

		private:

			double* scalar_p;

		};

	}
}

//distribution_function_template_D2Q9
namespace distribution_function_template_space {

	//D2Q9 distribution_functions_template--------------------------

	template <int X,int Y>
	class distribution_function_template_D2Q9  {
	public:

		distribution_function_template_D2Q9(double rho_initial = 0)
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
			// (i,j,q) <-- [i-1][j-1][q] <-- (i-1) * Y * 9 + (j - 1) * 9 + q
			int index = (i - 1) * Y * 9 + (j - 1) * 9 + q;
			return *(distribution_function_template_p + index);
		}

		//friend distribution_function_template_D2Q9<X, Y> operator+(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2);
		
		void blend(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2); //This function is used in the color gradient model to get the total distribution function.
		
		void streaming();

		//Solve for the equilibrium distribution function.
		//distribution_function_template_D2Q9.equilibrium(velocity2D_template,scalar_field);
		distribution_function_template_D2Q9<X,Y>& equilibrium(velocity2D_template_space::velocity2D_template<X,Y> &velocity,scalar_field_space::scalar_field<X,Y> &rho);

	protected:

		distribution_function_template_space::vector<double, 9> w;
		distribution_function_template_space::vector<double, 2> c[9];
		double* distribution_function_template_p;
		
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
	void distribution_function_template_D2Q9<X, Y>::blend(distribution_function_template_D2Q9<X,Y>& f1, distribution_function_template_D2Q9<X, Y>& f2) {
		for (int q = 0; q <= 8; q++)
			for (int i = 1; i <= X; i++)
				for (int j = 1; j <= Y; j++)
					(*this)(i, j, q) = f1(i, j, q) + f2(i, j, q);
		return;
	}

	template <int X, int Y>
	distribution_function_template_D2Q9<X,Y>& distribution_function_template_D2Q9<X, Y>::equilibrium( velocity2D_template_space::velocity2D_template<X, Y>& velocity,
		scalar_field_space::scalar_field<X, Y>& rho) {

		double uu = 0, cu = 0;
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				uu = velocity(i, j) * velocity(i, j);
				for (int q = 0; q <= 8; q++) {
					cu = velocity(i, j) * c[q];
					(*this)(i, j, q) = rho(i, j) * w(q) * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
				}
			}
		}
		return *this;
	}

}

#endif // !_DISTRIBUTION_FUNCTION_H_
