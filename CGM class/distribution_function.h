#ifndef _DISTRIBUTION_FUNCTION_H_
#define _DISTRIBUTION_FUNCTION_H_
#include <iostream>
#include <cmath>
#include <cstdio>
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
			if (i >= N || i < 0) { 
				printf("\nVector template is called, but the subscript exceeds the limit.\n"); 
			}
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
				if (i > X) {
					printf("\nvelocity2D_template is called, but the first subscript exceeds the limit: %d > X\n",i);
				}
				else if (i < 1) {
					printf("\nvelocity2D_template is called, but the first subscript exceeds the limit: %d < 1\n",i);
				}
				else if (j > Y) {
					printf("\nvelocity2D_template is called, but the second subscript exceeds the limit: %d > Y\n",j);
				}
				else if (j < 1) {
					printf("\nvelocity2D_template is called, but the second subscript exceeds the limit: %d < 1\n",j);
				}
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

	//state distribution_function_template_D2Q9 explicitly beforehand
	template<int X, int Y>
	class distribution_function_template_D2Q9;

	namespace scalar_field_space {

		//define a scalar_field template
		template<int X,int Y>
		class scalar_field {
		public:
			scalar_field(double scalar_initial = 0) : scalar_p(new double[X*Y])
			{
				for (int r = 0; r < X * Y; r++)
					*(scalar_p + r) = scalar_initial;
			}

			~scalar_field(){
				delete[]scalar_p;
			}

			virtual double& operator()(int i,int j) {
				// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
				if (i > X) {
					printf("\nscalar_field is called, but the first subscript exceeds the limit: %d > X\n", i);
				}
				else if (i < 1) {
					printf("\nscalar_field is called, but the first subscript exceeds the limit: %d < 1\n", i);
				}
				else if (j > Y) {
					printf("\nscalar_field is called, but the second subscript exceeds the limit: %d > Y\n", j);
				}
				else if (j < 1) {
					printf("\nscalar_field is called, but the second subscript exceeds the limit: %d < 1\n", j);
				}
				int index = (i - 1) * Y + (j - 1);
				return *(scalar_p + index);
			}

		protected:

			double* scalar_p;

		};

		//define a density_field template inheriting from scalar_field template
		template<int X,int Y>
		class density_field :public scalar_field<X, Y> {
		public:
			density_field(double density_initial = 1.0) :scalar_field<X, Y>(density_initial) {}

			virtual double& operator()(int i, int j) override{
				// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
				if (i > X) {
					printf("\ndensity_field is called, but the first subscript exceeds the limit: %d > X\n", i);
				}
				else if (i < 1) {
					printf("\ndensity_field is called, but the first subscript exceeds the limit: %d < 1\n", i);
				}
				else if (j > Y) {
					printf("\ndensity_field is called, but the second subscript exceeds the limit: %d > Y\n", j);
				}
				else if (j < 1) {
					printf("\ndensity_field is called, but the second subscript exceeds the limit: %d < 1\n", j);
				}
				int index = (i - 1) * Y + (j - 1);
				return *(scalar_field<X,Y>::scalar_p + index);
			}

			//density_field own function for calculate density from distribution function
			density_field& calculate(distribution_function_template_D2Q9<X, Y>& f) {
				double temp = 0;
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						temp = 0;
						for (int q = 0; q <= 8; q++) {
							temp += f(i, j, q);
						}
						(*this)(i, j) = temp;
					}
				}
				return *this;
			}
		};

		//define a phase_field template inheriting from scalar_field template
		template<int X, int Y>
		class phase_field :public scalar_field<X, Y> {
		public:
			phase_field(double phase_initial = 0) :scalar_field<X, Y>(phase_initial) {}
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
			if (i > X) {
				printf("\ndistribution_function_template_D2Q9 is called, but the first subscript exceeds the limit: %d > X\n", i);
			}
			else if (i < 1) {
				printf("\ndistribution_function_template_D2Q9 is called, but the first subscript exceeds the limit: %d < 1\n", i);
			}
			else if (j > Y) {
				printf("\ndistribution_function_template_D2Q9 is called, but the second subscript exceeds the limit: %d > Y\n", j);
			}
			else if (j < 1) {
				printf("\ndistribution_function_template_D2Q9 is called, but the second subscript exceeds the limit: %d < 1\n", j);
			}
			else if (q < 0) {
				printf("\ndistribution_function_template_D2Q9 is called, but the third subscript exceeds the limit: %d < 0\n", q);
			}
			else if (q > 8) {
				printf("\ndistribution_function_template_D2Q9 is called, but the third subscript exceeds the limit: %d > 8\n", q);
			}
			int index = (i - 1) * Y * 9 + (j - 1) * 9 + q;
			return *(distribution_function_template_p + index);
		}

		//friend distribution_function_template_D2Q9<X, Y> operator+(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2);
		
		void blend(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2); //This function is used in the color gradient model to get the total distribution function.
		
		void streaming();

		//Solve for the equilibrium distribution function.
		//distribution_function_template_D2Q9.equilibrium(velocity2D_template,scalar_field);
		distribution_function_template_D2Q9<X,Y>& equilibrium(velocity2D_template_space::velocity2D_template<X,Y> &velocity,scalar_field_space::scalar_field<X,Y> &rho);

		//Detect if the value of the distribution function is abnormal.
		bool detect();

	protected:

		distribution_function_template_space::vector<double, 9> w;
		distribution_function_template_space::vector<double, 2> c[9];
		double* distribution_function_template_p;
		
	};

	template <int X, int Y>
	inline void distribution_function_template_D2Q9<X, Y>::streaming() {
		double* temp_x = new double[X+1];
		double* temp_y = new double[Y+1];
		for (int q = 1; q <= 8; q++) {
			switch (q)
			{
			case 1:
				for (int j = 1; j <= Y; j++) temp_y[j] = this->operator()(X, j, q);
				for (int j = 1; j <= Y; j++) for (int i = X; i >= 2; i--) this->operator()(i, j, q) = this ->operator()(i - 1, j, q); 
				for (int j = 1; j <= Y; j++) this->operator()(1, j, q) = temp_y[j]; 
				break;
			case 2:
				for (int i = 1; i <= X; i++) temp_x[i] = this->operator()(i, Y, q);
				for (int i = 1; i <= X; i++) for (int j = Y; j >= 2; j--) this->operator()(i, j, q) = this ->operator()(i, j - 1, q); 
				for (int i = 1; i <= X; i++) this->operator()(i, 1, q) = temp_x[i]; 
				break;
			case 3:
				for (int j = 1; j <= Y; j++) temp_y[j] = this->operator()(1, j, q);
				for (int j = 1; j <= Y; j++) for (int i = 1; i <= X - 1; i++)  this->operator()(i, j, q) = this ->operator()(i + 1, j, q); 
				for (int j = 1; j <= Y; j++) this->operator()(X, j, q) = temp_y[j]; 
				break;
			case 4:
				for (int i = 1; i <= X; i++) temp_x[i] = this->operator()(i, 1, q);
				for (int i = 1; i <= X; i++) for (int j = 1; j <= Y - 1; j++) this->operator()(i, j, q) = this ->operator()(i, j + 1, q); 
				for (int i = 1; i <= X; i++) this->operator()(i, Y, q) = temp_x[i]; 
				break;
			case 5:
				for (int j = 1; j <= Y; j++) temp_y[j] = this->operator()(X, j, q);
				for (int i = 1; i <= X; i++) temp_x[i] = this->operator()(i, Y, q);
				for (int i = X; i >= 2; i--) for (int j = Y; j >= 2; j--) this->operator()(i, j, q) = this ->operator()(i - 1, j - 1, q); 
				for (int j = 2; j <= Y; j++) this->operator()(1, j, q) = temp_y[j - 1];
				for (int i = 2; i <= X; i++) this->operator()(i, 1, q) = temp_x[i - 1];
				this->operator()(1, 1, q) = temp_x[X];
				break;
			case 6:
				for (int i = 1; i <= X; i++) temp_x[i] = this->operator()(i, Y, q);
				for (int j = 1; j <= Y; j++) temp_y[j] = this->operator()(1, j, q);
				for (int i = 1; i <= X - 1; i++) for (int j = Y; j >= 2; j--) this->operator()(i, j, q) = this ->operator()(i + 1, j - 1, q); 
				for (int i = 1; i <= X - 1; i++) this->operator()(i, 1, q) = temp_x[i + 1];
				for (int j = 2; j <= Y; j++) this->operator()(X, j, q) = temp_y[j - 1];
				this->operator()(X, 1, q) = temp_x[1];
				break;
			case 7:
				for (int j = 1; j <= Y; j++) temp_y[j] = this->operator()(1, j, q);
				for (int i = 1; i <= X; i++) temp_x[i] = this->operator()(i, 1, q);
				for (int i = 1; i <= X - 1; i++) for (int j = 1; j <= Y - 1; j++) this->operator()(i, j, q) = this ->operator()(i + 1, j + 1, q); 
				for (int j = 1; j <= Y - 1; j++) this->operator()(X, j, q) = temp_y[j + 1];
				for (int i = 1; i <= X - 1; i++) this->operator()(i, Y, q) = temp_x[i + 1];
				this->operator()(X, Y, q) = temp_x[1];
				break;
			case 8:
				for (int j = 1; j <= Y; j++) temp_y[j] = this->operator()(X, j, q);
				for (int i = 1; i <= X; i++) temp_x[i] = this->operator()(i, 1, q);
				for (int i = X; i >= 2; i--) for (int j = 1; j <= Y - 1; j++) this->operator()(i, j, q) = this ->operator()(i - 1, j + 1, q); 
				for (int j = 1; j <= Y - 1; j++) this->operator()(1, j, q) = temp_y[j + 1];
				for (int i = 2; i <= X; i++) this->operator()(i, Y, q) = temp_x[i - 1];
				this->operator()(1, Y, q) = temp_x[X];
				break;
			default:
				break;
			}
		}
		delete[]temp_x;
		delete[]temp_y;
		return;
	}

	//template <int X,int Y>
	//distribution_function_template_D2Q9<X, Y> operator+(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2) {
	//}

	template <int X,int Y>
	inline void distribution_function_template_D2Q9<X, Y>::blend(distribution_function_template_D2Q9<X,Y>& f1, distribution_function_template_D2Q9<X, Y>& f2) {
		for (int q = 0; q <= 8; q++)
			for (int i = 1; i <= X; i++)
				for (int j = 1; j <= Y; j++)
					(*this)(i, j, q) = f1(i, j, q) + f2(i, j, q);
		return;
	}

	template <int X, int Y>
	inline distribution_function_template_D2Q9<X,Y>& distribution_function_template_D2Q9<X, Y>::equilibrium( velocity2D_template_space::velocity2D_template<X, Y>& velocity,
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

	template <int X, int Y>
	inline bool distribution_function_template_D2Q9<X, Y>::detect() {
		double min = 0, max = 0;
		int position[4]{ 1,1,1,1 };
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j < Y; j++) {
				for (int q = 0; q <= 8; q++) {
					if (isfinite((*this)(i, j, q)) != 0) {
						if ((*this)(i, j, q) >= 0) {
							min = this->operator()(i, j, q) < min ? (position[0] = i, position[1] = j, this->operator()(i, j, q)) : min;
							max = this->operator()(i, j, q) > max ? (position[2] = i, position[3] = j, this->operator()(i, j, q)) : max;
						}
						else {
							printf("\nThe value of distribution function (%d,%d,%d) is %.6f < 0.It doesn't satisfy the requirement that it must be positive.\n", i, j, q, (*this)(i, j, q));
							return false;
						}
					}
					else{ 
						printf("\nThe value of distribution function (%d,%d,%d) is infinite.\n", i, j, q);
						return false;
					}
				}
			}
		}
		printf("\nmax=%.6f at (%d,%d), min=%.6f at (%d,%d)\n", max, position[0], position[1], min, position[2], position[3]);
		return true;
	}

}

#endif // !_DISTRIBUTION_FUNCTION_H_
