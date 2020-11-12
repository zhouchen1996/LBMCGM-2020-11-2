#ifndef _DISTRIBUTION_FUNCTION_H_
#define _DISTRIBUTION_FUNCTION_H_

#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>
#include <string>

//vector template
namespace distribution_function_template_space {

	//important!!!!!!!!
	//vector template--------------------------------------
	template <typename T, int N>
	class vector {
	public:

		vector() :v(new T[N]){
			for (int i = 0; i < N; i++) {
				this->operator()(i) = T();
			}
		}

		vector(T initial) :v(new T[N])
		{
			for (int i = 0; i < N; i++) {
				this->operator()(i) = initial;
			}
		}

		vector(const vector<T, N>& a) : v(new T[N]) {
			for (int i = 0; i < N; i++)
				this->operator()(i) = a.v[i];
		}

		vector(const T(&a)[N]) : v(new T[N]) {
			for (int i = 0; i < N; i++)
				this->operator()(i) = a[i];
		}

		~vector() {
			delete[] v;
		}

		template <typename T, int N>
		friend T operator*(const vector<T, N>& a, const vector<T, N>& b);

		template <typename T, int N>
		friend vector<T, N> operator+(const vector<T, N>& a, const vector<T, N>& b);

		template <typename T, int N>
		friend vector<T, N> operator-(const vector<T, N>& a, const vector<T, N>& b);

		T& operator()(int i) {
			//start from 0
			if (i >= N || i < 0) { 
				printf("\nVector template is called, but the subscript exceeds the limit.\n"); 
			}
			return v[i];
		}

		T& operator()(int i) const {
			//start from 0
			if (i >= N || i < 0) {
				printf("\nVector template is called, but the subscript exceeds the limit.\n");
			}
			return v[i];
		}

		vector<T, N>& operator=(const vector<T, N>& a) {
			for (int i = 0; i < N; i++)
				//this->v[i] = a.v[i];
				this->operator()(i) = a(i);
			return *this;
		}

		vector<T, N>& operator=(const T(&a)[N]) {
			for (int i = 0; i < N; i++)
				this->operator()(i) = a[i];
			return *this;
		}

		vector<T, N>& operator+=(vector<T, N>& a) {
			for (int i = 0; i < N; i++)
				this->operator()(i) += this->operator()(i);
			return *this;
		}

	private:
		T* v;
	};

	// * for vector template
	template <typename T, int N>
	T operator*(const vector<T, N>& a, const vector<T, N>& b) {
		T temp(0);
		for (int i = 0; i < N; i++)
			temp += (a(i) * b(i));
		return temp;
	}

	// + for vector template 
	template <typename T, int N>
	vector<T, N> operator+(const vector<T, N>& a, const vector<T, N>& b) {
		vector<T, N> temp_vector;
		for (int i = 0; i < N; i++)
			temp_vector(i) = a(i) + b(i);
		return temp_vector;
	}

	// - for vector template 
	template <typename T, int N>
	vector<T, N> operator-(const vector<T, N>& a, const vector<T, N>& b) {
		vector<T, N> temp_vector;
		for (int i = 0; i < N; i++)
			temp_vector(i) = a(i) - b(i);
		return temp_vector;
	}

}

//matrix template
namespace distribution_function_template_space {

	//Different from vectors' indices starting from 0 , matrix indices start from 1.
	template <typename T, int X,int Y>
	class matrix {
	public:
		matrix() : matrix_p(new T[X * Y]){
			for (int r = 0; r < X * Y; r++)
				*(matrix_p + r) = T();
		}

		matrix(const T& initial) : matrix_p(new T[X * Y])
		{
			for (int r = 0; r < X * Y; r++)
				*(matrix_p + r) = initial;
		}

		matrix(const T(&a)[X][Y]) : matrix_p(new T[X * Y]) {
			for (int i = 1; i <= X; i++) {
				for (int j = 1; j <= Y; j++) {
					(*this)(i, j) = a[i - 1][j - 1];
				}
			}
		}

		matrix(const vector<T,X>&a) : matrix_p(new T[X * Y]) {
			//Construct a diagonal matrix
			for (int i = 1; i <= X; i++) {
				for (int j = 1; j <= Y; j++) {
					if (i == j) {
						(*this)(i, j) = a(i - 1);
					}
					else {
						(*this)(i, j) = 0;
					}
				}
			}
		}

		matrix(const matrix<T, X, Y>& matrix_) : matrix_p(new T[X * Y]) {
			for (int i = 1; i <= X; i++) {
				for (int j = 1; j <= Y; j++) {
					(*this)(i, j) = matrix_(i, j);
				}
			}
		}

		double& operator()(int i, int j) {
			// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
			if (i > X) {
				printf("\nmatrix is called, but the first subscript exceeds the limit: %d > X\n", i);
			}
			else if (i < 1) {
				printf("\nmatrix is called, but the first subscript exceeds the limit: %d < 1\n", i);
			}
			else if (j > Y) {
				printf("\nmatrix is called, but the second subscript exceeds the limit: %d > Y\n", j);
			}
			else if (j < 1) {
				printf("\nmatrix is called, but the second subscript exceeds the limit: %d < 1\n", j);
			}
			int index = (i - 1) * Y + (j - 1);
			return *(matrix_p + index);
		}

		double& operator()(int i, int j) const {
			// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
			if (i > X) {
				printf("\nmatrix is called, but the first subscript exceeds the limit: %d > X\n", i);
			}
			else if (i < 1) {
				printf("\nmatrix is called, but the first subscript exceeds the limit: %d < 1\n", i);
			}
			else if (j > Y) {
				printf("\nmatrix is called, but the second subscript exceeds the limit: %d > Y\n", j);
			}
			else if (j < 1) {
				printf("\nmatrix is called, but the second subscript exceeds the limit: %d < 1\n", j);
			}
			int index = (i - 1) * Y + (j - 1);
			return *(matrix_p + index);
		}

		matrix<T, X, Y>& operator=(const matrix<T, X, Y>& matrix_) {
			for (int i = 1; i <= X; i++) {
				for (int j = 1; j <= Y; j++) {
					(*this)(i, j) = matrix_(i, j);
				}
			}
			return *this;
		}

		matrix<T, X, Y>& operator+=(const matrix<T, X, Y>& matrix_) {
			for (int i = 1; i <= X; i++) {
				for (int j = 1; j <= Y; j++) {
					(*this)(i, j) += matrix_(i, j);
				}
			}
			return *this;
		}

		matrix<T, X, Y>& operator-=(const matrix<T, X, Y>& matrix_) {
			for (int i = 1; i <= X; i++) {
				for (int j = 1; j <= Y; j++) {
					(*this)(i, j) -= matrix_(i, j);
				}
			}
			return *this;
		}

		template<typename T, int X, int Y>
		friend matrix<T, X, Y> operator+(const matrix<T, X, Y>& matrix_1, const matrix<T, X, Y>& matrix_2);

		template<typename T, int X, int Y>
		friend matrix<T, X, Y> operator-(const matrix<T, X, Y>& matrix_1, const matrix<T, X, Y>& matrix_2);

		template<typename T, int X, int Y>
		friend vector<T, X> operator*(const matrix<T, X, Y>& matrix_, const vector<T, Y>& vector_);

		template<typename T, int X, int Y,int K>
		friend matrix<T, X, Y> operator*(const matrix<T, X, K>& matrix_1,const matrix<T, K, Y>& matrix_2);

		~matrix() {
			delete[]matrix_p;
		}

	protected:
		double* matrix_p;
	};

	template<typename T, int X, int Y>
	matrix<T, X, Y> operator+(const matrix<T, X, Y>& matrix_1, const matrix<T, X, Y>& matrix_2) {
		matrix<T, X, Y> matrix_temp(matrix_1);
		matrix_temp += matrix_2;
		return matrix_temp;
	}

	template<typename T, int X, int Y>
	matrix<T, X, Y> operator-(const matrix<T, X, Y>& matrix_1, const matrix<T, X, Y>& matrix_2) {
		matrix<T, X, Y> matrix_temp(matrix_1);
		matrix_temp -= matrix_2;
		return matrix_temp;
	}

	template<typename T, int X, int Y>
	vector<T, X> operator*(const matrix<T, X, Y>& matrix_, const vector<T, Y>& vector_) {
		vector<T, X> vector_temp;
		for (int i = 1; i <= X; i++) {
			vector_temp(i - 1) = 0;
			for (int j = 1; j <= Y; j++) {
				vector_temp(i - 1) += matrix_(i, j) * vector_(j - 1);
			}
		}
		return vector_temp;
	}

	template<typename T, int X, int Y, int K>
	matrix<T, X, Y> operator*(const matrix<T, X, K>& matrix_1, const matrix<T, K, Y>& matrix_2) {
		matrix<T, X, Y> matrix_temp;
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				matrix_temp(i, j) = 0;
				for (int k = 1; k <= K; k++) {
					matrix_temp(i, j) += matrix_1(i, k) * matrix_2(k, j);
				}
			}
		}
		return matrix_temp;
	}

}

//vector field 
namespace distribution_function_template_space {
	
	//state distribution_function_template_D2Q9 explicitly beforehand
	template<int X, int Y>
	class distribution_function_template_D2Q9;

	//define vector field
	template <int X, int Y>
	class vector2D_field {
	public:

		vector2D_field(double initial_vector_x = 0, double initial_vector_y = 0) :vector2D_field_p(new vector<double, 2>[X * Y]) {
			for (int r = 0; r < X * Y; r++)
				*(vector2D_field_p + r) = { initial_vector_x,initial_vector_y };
		}

		~vector2D_field() {
			delete[]vector2D_field_p;
		}

		vector<double, 2>& operator()(int i, int j) {
			// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
			if (i > X) {
				printf("\nvector2D_field is called, but the first subscript exceeds the limit: %d > X\n", i);
			}
			else if (i < 1) {
				printf("\nvector2D_field is called, but the first subscript exceeds the limit: %d < 1\n", i);
			}
			else if (j > Y) {
				printf("\nvector2D_field is called, but the second subscript exceeds the limit: %d > Y\n", j);
			}
			else if (j < 1) {
				printf("\nvector2D_field is called, but the second subscript exceeds the limit: %d < 1\n", j);
			}
			int index = (i - 1) * Y + (j - 1);
			return *(vector2D_field_p + index);
		}

	protected:
		vector<double, 2>* vector2D_field_p;
	};


	namespace vector2D_field_space {

		template<int X,int Y>
		class force2D_field :public vector2D_field<X, Y> {
		public:
			force2D_field(double initial_force2D_x = 0, double initial_force2D_y = 0) :vector2D_field<X, Y>(initial_force2D_x, initial_force2D_y) {}

			vector<double, 2>& operator()(int i, int j) {
				// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
				if (i > X) {
					printf("\nforce2D_field is called, but the first subscript exceeds the limit: %d > X\n", i);
				}
				else if (i < 1) {
					printf("\nforce2D_field is called, but the first subscript exceeds the limit: %d < 1\n", i);
				}
				else if (j > Y) {
					printf("\nforce2D_field is called, but the second subscript exceeds the limit: %d > Y\n", j);
				}
				else if (j < 1) {
					printf("\nforce2D_field is called, but the second subscript exceeds the limit: %d < 1\n", j);
				}
				int index = (i - 1) * Y + (j - 1);
				return *(this->vector2D_field_p + index);
			}
		};

	}

	//state scalar_field and density_field explicitly beforehand
	namespace  scalar_field_space {

		template<int X, int Y>
		class scalar_field;

		template<int X,int Y>
		class density_field;
	}

	namespace vector2D_field_space {

		template <int X, int Y>
		class velocity2D_field :public vector2D_field<X,Y>{
		public:

			velocity2D_field(double initial_velocity_x = 0, double initial_velocity_y = 0) :vector2D_field<X, Y>(initial_velocity_x, initial_velocity_y){}

			vector<double, 2>& operator()(int i, int j) {
				// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
				if (i > X) {
					printf("\nvelocity2D_field is called, but the first subscript exceeds the limit: %d > X\n",i);
				}
				else if (i < 1) {
					printf("\nvelocity2D_field is called, but the first subscript exceeds the limit: %d < 1\n",i);
				}
				else if (j > Y) {
					printf("\nvelocity2D_field is called, but the second subscript exceeds the limit: %d > Y\n",j);
				}
				else if (j < 1) {
					printf("\nvelocity2D_field is called, but the second subscript exceeds the limit: %d < 1\n",j);
				}
				int index = (i - 1) * Y + (j - 1);
				return *(this->vector2D_field_p + index);
			}
			
			//velocity2D_field own function for calculate velocity from distribution function
			velocity2D_field<X,Y>& calculate(distribution_function_template_D2Q9<X, Y>& f, scalar_field_space::density_field<X,Y>& rho) {
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						this->operator()(i, j)(0) = (f(i, j, 1) + f(i, j, 5) + f(i, j, 8) - f(i, j, 3) - f(i, j, 6) - f(i, j, 7))/rho(i,j);
						this->operator()(i, j)(1) = (f(i, j, 2) + f(i, j, 5) + f(i, j, 6) - f(i, j, 4) - f(i, j, 7) - f(i, j, 8))/rho(i,j);
					}
				}
				return *this;
			}

			//(Force term Guo.et al)velocity2D_field own function for calculate velocity from distribution function
			velocity2D_field<X, Y>& calculate(distribution_function_template_D2Q9<X, Y>& f, scalar_field_space::density_field<X, Y>& rho,vector2D_field_space::force2D_field<X,Y>& force2D) {
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						this->operator()(i, j)(0) = (f(i, j, 1) + f(i, j, 5) + f(i, j, 8) - f(i, j, 3) - f(i, j, 6) - f(i, j, 7) + 0.5 * force2D(i, j)(0)) / rho(i, j);
						this->operator()(i, j)(1) = (f(i, j, 2) + f(i, j, 5) + f(i, j, 6) - f(i, j, 4) - f(i, j, 7) - f(i, j, 8) + 0.5 * force2D(i, j)(1)) / rho(i, j);
					}
				}
				return *this;
			}

			//velocity2D_field own function for judge if the velocity value is normal
			bool detect() {
				double min = 0, max = 0;
				double temp = 0;
				int position[4]{ 1,1,1,1 };
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j < Y; j++) {
						if ((isfinite((*this)(i, j)(0)) && isfinite((*this)(i, j)(1))) != 0) {
							temp = sqrt((*this)(i, j)(0) * (*this)(i, j)(0) + (*this)(i, j)(1) * (*this)(i, j)(1));
							min = temp < min ? (position[0] = i, position[1] = j, temp) : min;
							max = temp > max ? (position[2] = i, position[3] = j, temp) : max;
						}
						else {
							printf("\nvelocity at (%d,%d) is infinite.\n", i, j);
							return false;
						}
					}
				}
				printf("\nmax=%.6f at (%d,%d), min=%.6f at (%d,%d)\n", max, position[0], position[1], min, position[2], position[3]);
				return true;
			}

		};

	}

}

//scalar_field  template
namespace distribution_function_template_space {


	//state distribution_function_template_D2Q9 explicitly beforehand
	template<int X, int Y>
	class distribution_function_template_D2Q9;

	//state areatype and fluidtype explicitly beforehand
	enum class areatype;

	enum class fluidtype;

	//state area_field and fluid_field explicitly beforehand
	template<typename T, int X, int Y>
	class area_field;

	template<typename T, int X, int Y>
	class fluid_field;

	namespace scalar_field_space {

		//define a scalar_field template
		template<int X,int Y>
		class scalar_field {
		public:
			scalar_field(double scalar_initial = 0) : scalar_p(new double[(X + 2) * (Y + 2)])
			{
				for (int r = 0; r < (X + 2) * (Y + 2); r++)
					*(scalar_p + r) = scalar_initial;
			}

			~scalar_field(){
				delete[]scalar_p;
			}

			double& operator()(int i,int j) {
				// (i,j) <-- [i][j] <-- i*(Y+2) + j
				if (i > X+1) {
					printf("\nscalar_field is called, but the first subscript exceeds the limit: %d > X+1\n", i);
				}
				else if (i < 0) {
					printf("\nscalar_field is called, but the first subscript exceeds the limit: %d < 0\n", i);
				}
				else if (j > Y+1) {
					printf("\nscalar_field is called, but the second subscript exceeds the limit: %d > Y+1\n", j);
				}
				else if (j < 0) {
					printf("\nscalar_field is called, but the second subscript exceeds the limit: %d < 0\n", j);
				}
				int index = i * (Y + 2) + j;
				return *(scalar_p + index);
			}

			//solve the sum of two scalar fields
			scalar_field<X, Y>& blend(scalar_field<X,Y>&scalar_1, scalar_field<X, Y>& scalar_2) {
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						(*this)(i, j) = scalar_1(i, j) + scalar_2(1, 2);
					}
				}
				return *this;
			}

		protected:

			double* scalar_p;

		};

		//define a density_field template inheriting from scalar_field template
		template<int X,int Y>
		class density_field :public scalar_field<X, Y> {
		public:

			density_field(double density_initial = 1.0, double solid_density_ = 1.0) :scalar_field<X, Y>(density_initial), solid_density(solid_density_) {}
			
			//set virtual solid density
			density_field(area_field<areatype,X,Y>& areafield, double density_initial = 1.0, double solid_density_ = 1.0) :scalar_field<X, Y>(density_initial), solid_density(solid_density_) {
				
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						if (areafield(i, j) == areatype::SL || areafield(i, j) == areatype::SB) {
							(*this)(i, j) = solid_density;
						}
					}
				}

			}

			double& operator()(int i, int j) {
				// (i,j) <-- [i][j] <-- i*(Y+2) + j
				if (i > X + 1) {
					printf("\ndensity_field is called, but the first subscript exceeds the limit: %d > X+1\n", i);
				}
				else if (i < 0) {
					printf("\ndensity_field is called, but the first subscript exceeds the limit: %d < 0\n", i);
				}
				else if (j > Y + 1) {
					printf("\ndensity_field is called, but the second subscript exceeds the limit: %d > Y+1\n", j);
				}
				else if (j < 0) {
					printf("\ndensity_field is called, but the second subscript exceeds the limit: %d < 0\n", j);
				}
				int index = i * (Y + 2) + j;
				return *(this->scalar_p + index);
			}

			//carefully! 谨慎使用
			//for calculating phase field gradient
			//When the index may exceeds the limit,the density is virtual solid density(a functor convenient for setting wettability conditions).
			double operator()(int i, int j, int) {
				// (i,j) <-- [i][j] <-- i*(Y+2) + j
				if (i > X || i < 1 || j > Y || j < 1) {
					return solid_density;
				}
				int index = i * (Y + 2) + j;
				return *(this->scalar_p + index);
			}

			density_field<X, Y>& set_virtual_solid_density(double rho_s) {
				solid_density = rho_s;
				return *this;
			}

			//density_field own function for calculate density from distribution function
			density_field<X,Y>& calculate(distribution_function_template_D2Q9<X, Y>& f) {
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

			//density_field own function for judge if the density value is normal
			bool detect() {
				double min = 0, max = 0;
				int position[4]{ 1,1,1,1 };
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j < Y; j++) {
						if (isfinite((*this)(i, j)) != 0) {
							if ((*this)(i, j) >= 0) {
								min = this->operator()(i, j) < min ? (position[0] = i, position[1] = j, this->operator()(i, j)) : min;
								max = this->operator()(i, j) > max ? (position[2] = i, position[3] = j, this->operator()(i, j)) : max;
							}
							else {
								printf("\ndensity at (%d,%d) is %.6f < 0.It doesn't satisfy the requirement that it must be positive.\n", i, j, (*this)(i, j));
								return false;
							}
						}
						else {
							printf("\ndensity at (%d,%d) is infinite.\n", i, j);
							return false;
						}
					}
				}
				printf("\nmax=%.6f at (%d,%d), min=%.6f at (%d,%d)\n", max, position[0], position[1], min, position[2], position[3]);
				return true;
			}

		protected:

			double solid_density;

		};

		//define a phase_field template inheriting from scalar_field template
		template<int X, int Y>
		class phase_field :public scalar_field<X, Y> {
		public:

			phase_field(double phase_initial = 0,double solid_phase_field_ = 0.0) :scalar_field<X, Y>(phase_initial),solid_phase_field(solid_phase_field_) {}

			double& operator()(int i, int j) {
				// (i,j) <-- [i][j] <-- i*(Y+2) + j
				if (i > X+1) {
					printf("\nphase_field is called, but the first subscript exceeds the limit: %d > X+1\n", i);
				}
				else if (i < 0) {
					printf("\nphase_field is called, but the first subscript exceeds the limit: %d < 0\n", i);
				}
				else if (j > Y+1) {
					printf("\nphase_field is called, but the second subscript exceeds the limit: %d > Y+1\n", j);
				}
				else if (j < 0) {
					printf("\nphase_field is called, but the second subscript exceeds the limit: %d < 0\n", j);
				}
				int index = i * (Y + 2) + j;
				return *(this->scalar_p + index);
			}

			//carefully! 谨慎使用
			//for calculating phase field gradient
			//When calculating the gradient, the index may exceeds the limit(a functor convenient for setting wettability conditions).
			double operator()(int i, int j,int) {
				// (i,j) <-- [i][j] <-- i*(Y+2) + j
				if (i > X || i < 1 || j > Y || j < 1) {
					return solid_phase_field;
				}
				int index = i * (Y + 2) + j;
				return *(this->scalar_p + index);
			}

			//These member function is for the color gradient model.
			//--1--phase_field's own function to calculate phase field (r1-r2)/(r1+r2) or (r1/rr1-r2/rr2)/(r1/rr1+r2/rr2)
			phase_field<X, Y>& calculate(scalar_field_space::density_field<X,Y>& rho_1, scalar_field_space::density_field<X, Y>& rho_2,double rho_reference_1 = 1, double rho_refernece_2 = 1) {
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						double a = rho_1(i, j) / rho_reference_1 - rho_2(i, j) / rho_refernece_2;
						double b = rho_1(i, j) / rho_reference_1 + rho_2(i, j) / rho_refernece_2;
						if ( b > 0) {
							(*this)(i, j) = a / b;
						}
						else {
							(*this)(i, j) = a / b;
							printf("--1--phase_field<X, Y>::calculate(rho1,rho2,1,1) is called,but at (%d,%d) rho1 + rho2 = %.5f< 0", i, j, b);
						}
					}
				}
				return *this;
			}

			//--2--phase_field's own function to calculate phase field and limit the scope to a certain fluid region F,FB,FL,inlet,outlet,inlet_S,outlet_S except S,SB,SL 
			//(r1-r2)/(r1+r2) or (r1/rr1-r2/rr2)/(r1/rr1+r2/rr2)
			phase_field<X, Y>& calculate(scalar_field_space::density_field<X, Y>& rho_1, scalar_field_space::density_field<X, Y>& rho_2, area_field<areatype,X,Y>& areafield,double rho_reference_1 = 1, double rho_refernece_2 = 1) {
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						if (areafield(i, j) != areatype::SL && areafield(i, j) != areatype::SB) {
							double a = rho_1(i, j) / rho_reference_1 - rho_2(i, j) / rho_refernece_2;
							double b = rho_1(i, j) / rho_reference_1 + rho_2(i, j) / rho_refernece_2;
							if (b > 0) {
								(*this)(i, j) = a / b;
							}
							else {
								(*this)(i, j) = a / b;
								printf("--2--phase_field<X, Y>::calculate(rho1,rho2,areafield,1,1) is called,but at (%d,%d) rho1 + rho2 = %.5f< 0", i, j, b);
							}
						}
					}
				}
				return *this;
			}

			//--3--phase_field's own function to calculate phase field (r1-r2)/r
			phase_field<X, Y>& calculate(scalar_field_space::density_field<X, Y>& rho_1, scalar_field_space::density_field<X, Y>& rho_2, scalar_field_space::density_field<X, Y>& rho) {
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						(*this)(i, j) = (rho_1(i, j) - rho_2(i, j)) / rho(i,j);
					}
				}
				return *this;
			}

			//--1--distinguish the fluidfield by the phase field and delta
			phase_field<X, Y>& find_interface_by_phasefield(fluid_field<fluidtype,X,Y>& fluidfield, area_field<areatype,X,Y>& areafield,double delta = 0.7){
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						if (areafield(i, j) != areatype::SL && areafield(i, j) != areatype::SB) {
							if (isfinite((*this)(i, j)) != 0) {
								if ((*this)(i, j) > delta) {
									fluidfield = fluidtype::red;
								}
								else if ((*this)(i, j) < -delta) {
									fluidfield = fluidtype::blue;
								}
								else
								{
									fluidfield = fluidtype::interface;
								}
							}
							else {
								printf("--1--phase_field<X, Y>::find_interface_by_phasefield() is called,but ths value of phase_field at (%d,%d) is infinite.", i, j);
								return *this;
							}
						}
					}
				}
				return *this;
			}

			//--2--distinguish the fluidfield by the phase field GRADIENT and delta_2
			phase_field<X, Y>& find_interface_by_gradient(fluid_field<fluidtype, X, Y>& fluidfield, area_field<areatype,X,Y>& areafield,double delta_2 = 0.05) {
				
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						if (areafield(i, j) != areatype::SL && areafield(i, j) != areatype::SB) {
							if (isfinite((*this)(i, j)) != 0) {
								if (phase_field_gradient(i, j) * phase_field_gradient(i, j) > delta_2) {
									fluidfield = fluidtype::interface;
								}
								else if ((*this)(i, j) >0) {
									fluidfield = fluidtype::red;
								}
								else
								{
									fluidfield = fluidtype::blue;
								}
							}
							else {
								printf("--2--phase_field<X, Y>::find_interface_by_gradient() is called,but ths value of phase_field at (%d,%d) is infinite.", i, j);
								return *this;
							}
						}
					}
				}

				return *this;
			}
			
			//this is traditional wetting boundary condtion ,set virtual solid density or virtual solid phase field
			//solve for the phase field gradient( use the special functor (int i, int j, int=0) )
			//--1--
			phase_field<X, Y>& calculate_phase_field_gradient(area_field<areatype, X, Y>& areafield) {
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						if (areafield(i, j) == areatype::FL || areafield(i, j) == areatype::FB) {

							//phase_field_gradient(i, j)(0) = 1.0 / 12.0 * ((*this)(i + 1, j + 1) - (*this)(i - 1, j + 1))
							//	+ 1.0 / 3.0 * ((*this)(i + 1, j) - (*this)(i - 1, j)) + 1.0 / 12.0 * ((*this)(i + 1, j - 1) - (*this)(i - 1, j - 1));
							//phase_field_gradient(i, j)(1) = 1.0 / 12.0 * ((*this)(i + 1, j + 1) - (*this)(i + 1, j - 1))
							//	+ 1.0 / 3.0 * ((*this)(i, j + 1) - (*this)(i, j - 1)) + 1.0 / 12.0 * ((*this)(i - 1, j + 1) - (*this)(i - 1, j - 1));
							
							phase_field_gradient(i, j)(0) = 1.0 / 12.0 * ((*this)(i + 1, j + 1, 0) - (*this)(i - 1, j + 1, 0))
								+ 1.0 / 3.0 * ((*this)(i + 1, j, 0) - (*this)(i - 1, j, 0)) + 1.0 / 12.0 * ((*this)(i + 1, j - 1, 0) - (*this)(i - 1, j - 1, 0));
							
							phase_field_gradient(i, j)(1) = 1.0 / 12.0 * ((*this)(i + 1, j + 1, 0) - (*this)(i + 1, j - 1, 0))
								+ 1.0 / 3.0 * ((*this)(i, j + 1, 0) - (*this)(i, j - 1, 0)) + 1.0 / 12.0 * ((*this)(i - 1, j + 1, 0) - (*this)(i - 1, j - 1, 0));

						}
					}
				}
				return *this;
			}

			//this is new wetting boundary condtion ,first set phase field gradient at FL,then set FB's phase field gradient
			//solve for the phase field gradient
			//--2--
			phase_field<X, Y>& calculate_phase_field_gradient_new(area_field<areatype, X, Y>& areafield) {
				for (int i = 2; i <= X-1; i++) {
					for (int j = 2; j <= Y-1; j++) {
						if (areafield(i, j) == areatype::FL) {

							phase_field_gradient(i, j)(0) = 1.0 / 12.0 * ((*this)(i + 1, j + 1) - (*this)(i - 1, j + 1))
								+ 1.0 / 3.0 * ((*this)(i + 1, j) - (*this)(i - 1, j)) + 1.0 / 12.0 * ((*this)(i + 1, j - 1) - (*this)(i - 1, j - 1));
							phase_field_gradient(i, j)(1) = 1.0 / 12.0 * ((*this)(i + 1, j + 1) - (*this)(i + 1, j - 1))
								+ 1.0 / 3.0 * ((*this)(i, j + 1) - (*this)(i, j - 1)) + 1.0 / 12.0 * ((*this)(i - 1, j + 1) - (*this)(i - 1, j - 1));

						}
					}
				}
				return *this;
			}

			//visit phase_field_gradient
			vector<double, 2>& gradient(int i, int j) {
				return phase_field_gradient(i, j);
			}

			protected:
				double solid_phase_field;
				vector2D_field<X, Y> phase_field_gradient;
		};

	}

}

//define classes associated with area type and fluid type that represents different areas
namespace distribution_function_template_space {

	enum class areatype{S,F,SB,SL,FB,FL,inlet,INLET,outlet,OUTLET,FB_inlet,FB_outlet};
	enum class fluidtype{red,blue,interface};

	template<typename T,int X, int Y>
	class area_field {
	public:

		area_field() :area_field_p(new T[(X + 2) * (Y + 2)]) {}

		area_field(T area_initial_type) :area_field_p(new T[(X + 2) * (Y + 2)]) {
			for (int r = 0; r < (X + 2) * (Y + 2); r++) {
				*(area_field_p + r) = area_initial_type;
			}
		}
		~area_field() {
			delete[]area_field_p;
		}
		T& operator()(int i,int j) {
			// (i,j) <-- [i][j] <-- i * (Y + 2) + j
			if (i > X + 1) {
				printf("\narea_field is called, but the first subscript exceeds the limit: %d > X + 1\n", i);
			}
			else if (i < 0) {
				printf("\narea_field is called, but the first subscript exceeds the limit: %d < 0\n", i);
			}
			else if (j > Y + 1) {
				printf("\narea_field is called, but the second subscript exceeds the limit: %d > Y + 1\n", j);
			}
			else if (j < 0) {
				printf("\narea_field is called, but the second subscript exceeds the limit: %d < 0\n", j);
			}
			int index = i * (Y + 2) + j;
			return *(area_field_p + index);
		}
	protected:
		T* area_field_p;
	};

	template<typename T, int X, int Y>
	class fluid_field {
	public:

		fluid_field() :fluid_field_p(new T[X * Y]) {}

		fluid_field(T fluid_initial_type) :fluid_field_p(new T[X * Y]) {
			for (int r = 0; r < X * Y; r++) {
				*(fluid_field_p + r) = fluid_initial_type;
			}
		}
		~fluid_field() {
			delete[]fluid_field_p;
		}
		T& operator()(int i, int j) {
			// (i,j) <-- [i-1][j-1] <-- (i-1) * Y + (j-1)
			if (i > X) {
				printf("\nfluid_field is called, but the first subscript exceeds the limit: %d > X\n", i);
			}
			else if (i < 1) {
				printf("\nfluid_field is called, but the first subscript exceeds the limit: %d < 1\n", i);
			}
			else if (j > Y) {
				printf("\nfluid_field is called, but the second subscript exceeds the limit: %d > Y\n", j);
			}
			else if (j < 1) {
				printf("\nfluid_field is called, but the second subscript exceeds the limit: %d < 1\n", j);
			}
			int index = (i - 1) * Y + (j - 1);
			return *(fluid_field_p + index);
		}
	protected:
		T* fluid_field_p;
	};

}

//distribution_function_template_D2Q9
namespace distribution_function_template_space {

	//D2Q9 distribution_functions_template--------------------------

	template <int X, int Y>
	class distribution_function_template_D2Q9 {
	public:

		distribution_function_template_D2Q9(double rho_initial = 1.0, double nu_ = 1.0 / 6.0)
			: distribution_function_template_p(new double[X * Y * 9]), density(rho_initial), nu(nu_), S({ 1.0,1.64,1.54,1.0,1.2,1.0,1.2, 1.0 / (3.0 * nu_ + 0.5), 1.0 / (3.0 * nu_ + 0.5) })
		{
			for (int i = 1; i <= X; i++)
				for (int j = 1; j <= Y; j++)
					for (int q = 0; q <= 8; q++)
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

		//use carefully!
		//for boundary condition
		double& operator()(int i, int j, int q,int) {
			// (i,j,q) <-- [i-1][j-1][q] <-- (i-1) * Y * 9 + (j - 1) * 9 + q
			if (i == X + 1) {
				i = 1;
			}
			if (i == 0) {
				i = X;
			}
			if (j == Y + 1) {
				j = 1;
			}
			if (j == 0) {
				j = Y;
			}
			if (q < 0) {
				printf("\ndistribution_function_template_D2Q9 is called, but the third subscript exceeds the limit: %d < 0\n", q);
			}
			else if (q > 8) {
				printf("\ndistribution_function_template_D2Q9 is called, but the third subscript exceeds the limit: %d > 8\n", q);
			}
			int index = (i - 1) * Y * 9 + (j - 1) * 9 + q;
			return *(distribution_function_template_p + index);
		}

		
		distribution_function_template_D2Q9<X, Y>& streaming();

		//Solve for the equilibrium distribution function.
		//distribution_function_template_D2Q9.equilibrium(velocity2D_field,scalar_field);
		distribution_function_template_D2Q9<X, Y>& equilibrium(vector2D_field_space::velocity2D_field<X, Y>& velocity, scalar_field_space::scalar_field<X, Y>& rho);
		distribution_function_template_D2Q9<X, Y>& equilibrium();
		double equilibrium(int i, int j, int q);

		//Detect if the value of the distribution function is abnormal.
		bool detect();

		//--1--SRT
		distribution_function_template_D2Q9<X, Y>& single_phase_collison_SRT();
		//--2--MRT
		distribution_function_template_D2Q9<X, Y>& single_phase_collison_MRT();

		//boundary condition : bounce back
		distribution_function_template_D2Q9<X, Y>& bounce_back_halfway();

		//boundary condition : velocity inlet
		distribution_function_template_D2Q9<X, Y>& inlet_velocity(vector<double, 2> u_inlet);

		//boundary condition : pressure inlet
		//distribution_function_template_D2Q9<X, Y>& inlet_pressure(vector<double, 2> u_inlet);

		//boundary condition : outlflow oulet (Neumann boundary condition)
		//"Evaluation of outflow boundary conditions for two-phase lattice Boltzmann equation"(2013)
		distribution_function_template_D2Q9<X, Y>& oulet_outlflow_Neumann();

		//boundary condition : outlflow oulet (convective boundary condition)
		//"Evaluation of outflow boundary conditions for two-phase lattice Boltzmann equation"(2013)
		distribution_function_template_D2Q9<X, Y>& oulet_outlflow_convective();

	protected:

		static vector<double, 9> w;
		static vector<double, 2> c[9];
		double* distribution_function_template_p;

	public:

		//It works for all distribution functions。
		static vector2D_field_space::velocity2D_field<X, Y> velocity;
		static area_field<areatype, X, Y> area;
		static matrix<double, 9, 9> M;
		static matrix<double, 9, 9> InM;

		static void setarea(std::string& area_file_name);

		//It varies depending on the different distribution functions.
		double nu;
		vector<double, 9> S;
		scalar_field_space::density_field<X, Y> density;
	};

	//initialize (constant)
	//--1--
	template <int X, int Y>
	vector<double, 9> distribution_function_template_D2Q9<X, Y>::w({ 4.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0 });

	//--2--
	template <int X, int Y>
	vector<double, 2> distribution_function_template_D2Q9<X, Y>::c[9]
	{ vector<double, 2>({0, 0}),vector<double, 2>({1, 0}),vector<double, 2>({0, 1}),
		vector<double, 2>({-1, 0}),vector<double, 2>({0, -1}),vector<double, 2>({1, 1}),
		vector<double, 2>({-1, 1}),vector<double, 2>({-1, -1}),vector<double, 2>({1, -1}) };

	//--3--
	template <int X, int Y>
	matrix<double, 9, 9> distribution_function_template_D2Q9<X, Y>::M({
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1},
	{-4,-1,-1,-1,-1, 2, 2, 2, 2},
	{ 4,-2,-2,-2,-2, 1, 1, 1, 1},
	{ 0, 1, 0,-1, 0, 1,-1,-1, 1},
	{ 0,-2, 0, 2, 0, 1,-1,-1, 1},
	{ 0, 0, 1, 0,-1, 1, 1,-1,-1},
	{ 0, 0,-2, 0, 2, 1, 1,-1,-1},
	{ 0, 1,-1, 1,-1, 0, 0, 0, 0},
	{ 0, 0, 0, 0, 0, 1,-1, 1,-1} });

	//--4--
	template <int X, int Y>
	matrix<double, 9, 9> distribution_function_template_D2Q9<X, Y>::InM({
	{1.0 / 9.0, -1.0 / 9.0,  1.0 / 9.0,    0.0,    0.0,    0.0,   0.0,   0.0,  0.0},
	{1.0 / 9.0, -1.0 / 36.0,-1.0 / 18.0,1.0 / 6.0,-1.0 / 6.0,0.0,0.0,1.0 / 4.0,0.0},
	{1.0 / 9.0, -1.0 / 36.0,-1.0 / 18.0,0.0,0.0,1.0 / 6.0,-1.0 / 6.0,-1.0 / 4.0,0.0},
	{1.0 / 9.0, -1.0 / 36.0,-1.0 / 18.0,-1.0 / 6.0,1.0 / 6.0,0.0,0.0,1.0 / 4.0,0.0},
	{1.0 / 9.0, -1.0 / 36.0,-1.0 / 18.0,0.0,0.0,-1.0 / 6.0,1.0 / 6.0,-1.0 / 4.0,0.0},
	{1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0,1.0 / 6.0,1.0 / 12.0,1.0 / 6.0,1.0 / 12.0,0.0,1.0 / 4.0},
	{1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0,-1.0 / 6.0,-1.0 / 12.0,1.0 / 6.0,1.0 / 12.0,0.0,-1.0 / 4.0},
	{1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0,-1.0 / 6.0,-1.0 / 12.0,-1.0 / 6.0,-1.0 / 12.0,0.0,1.0 / 4.0},
	{1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0,1.0 / 6.0,1.0 / 12.0,-1.0 / 6.0,-1.0 / 12.0,0.0,-1.0 / 4.0} });

	//initialize (optional)
	//--1--
	template <int X, int Y>
	vector2D_field_space::velocity2D_field<X, Y> distribution_function_template_D2Q9<X, Y>::velocity(0.001, 0);

	//--2--
	template<int X, int Y>
	area_field<areatype, X, Y> distribution_function_template_D2Q9<X, Y>::area(areatype::F);

	template <int X, int Y>
	inline distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::streaming() {
		double* temp_x = new double[X + 1];
		double* temp_y = new double[Y + 1];
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
		return *this;
	}


	//several equilibrium functions
	//--1--
	template <int X, int Y>
	inline distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::equilibrium(vector2D_field_space::velocity2D_field<X, Y>& velocity,
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
	//--2--
	template <int X, int Y>
	inline distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::equilibrium() {

		double uu = 0, cu = 0;
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				uu = velocity(i, j) * velocity(i, j);
				for (int q = 0; q <= 8; q++) {
					cu = velocity(i, j) * c[q];
					(*this)(i, j, q) = density(i, j) * w(q) * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
				}
			}
		}

		return *this;
	}
	//--3--
	template <int X, int Y>
	inline double distribution_function_template_D2Q9<X, Y>::equilibrium(int i, int j, int q) {
		double uu = velocity(i, j) * velocity(i, j), cu = velocity(i, j) * c[q];
		return density(i, j) * w(q) * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
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
					else {
						printf("\nThe value of distribution function (%d,%d,%d) is infinite.\n", i, j, q);
						return false;
					}
				}
			}
		}
		printf("\nmax=%.6f at (%d,%d), min=%.6f at (%d,%d)\n", max, position[0], position[1], min, position[2], position[3]);
		return true;
	}

	//SRT
	template <int X, int Y>
	distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::single_phase_collison_SRT() {
		double omega = 1.0 / (3.0 * nu + 0.5);
		for (int i = 1; i <= X; i++)
			for (int j = 1; j <= Y; j++)
				for (int q = 0; q <= 8; q++)
					(*this)(i, j, q) = (1 - omega) * (*this)(i, j, q) + omega * (*this).equilibrium(i, j, q);
		return *this;
	}

	//MRT
	template <int X, int Y>
	distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::single_phase_collison_MRT() {
		vector<double, 9> moment_pre_collison, moment_eq, moment_post_collison, f_vector_pre_collison, f_vector_eq, f_vector_post_collison;
		matrix<double, 9, 9> I(vector<double, 9>(1.0)), Matrix_S(S);
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {

				for (int q = 0; q <= 8; q++)
					f_vector_pre_collison(q) = (*this)(i, j, q);
				
				for (int q = 0; q <= 8; q++)
					f_vector_eq(q) = (*this).equilibrium(i, j, q);

				//Convert to the moment space
				moment_pre_collison = M * f_vector_pre_collison;

				//Calculate the temporary equilibrium moment
				moment_eq = M * f_vector_eq;

				//Collision in moment space
				moment_post_collison = (I - Matrix_S) * moment_pre_collison + Matrix_S * moment_eq;

				//Back to particle space
				f_vector_post_collison = InM * moment_post_collison;

				for (int q = 0; q <= 8; q++)
					(*this)(i, j, q) = f_vector_post_collison(q);
			}
		}		
		return *this;
	}

	template <int X, int Y>
	distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::bounce_back_halfway() {
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (area(i,j) == areatype::FB || area(i, j) == areatype::FB_inlet) {
					
					if (area(i - 1,j) == areatype::SB) {
						(*this)(i, j, 1, 0) = (*this)(i - 1, j, 3, 0);
					}
					if (area(i, j - 1) == areatype::SB) {
						(*this)(i, j, 2, 0) = (*this)(i, j - 1, 4, 0);
					}
					if (area(i + 1, j) == areatype::SB) {
						(*this)(i, j, 3, 0) = (*this)(i + 1, j, 1, 0);
					}
					if (area(i, j + 1) == areatype::SB) {
						(*this)(i, j, 4, 0) = (*this)(i, j + 1, 2, 0);
					}
					if (area(i - 1, j - 1) == areatype::SB) {
						(*this)(i, j, 5, 0) = (*this)(i - 1, j - 1, 7, 0);
					}
					if (area(i + 1, j - 1) == areatype::SB) {
						(*this)(i, j, 6, 0) = (*this)(i + 1, j - 1, 8, 0);
					}
					if (area(i + 1, j + 1) == areatype::SB) {
						(*this)(i, j, 7, 0) = (*this)(i + 1, j + 1, 5, 0);
					}
					if (area(i - 1, j + 1) == areatype::SB) {
						(*this)(i, j, 8, 0) = (*this)(i - 1, j + 1, 6, 0);
					}

				}
			}
		}
		return *this;
	}

	template <int X, int Y>
	distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::inlet_velocity(vector<double,2> u_inlet) {
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (area(i, j) == areatype::inlet || area(i, j) == areatype::FB_inlet) {

					if (area(i - 1, j) == areatype::INLET) {
						(*this)(i, j, 1, 0) = (*this)(i - 1, j, 3, 0) - 2.0 / 3.0 * (*this).density(i, j) * (c[3] * u_inlet);
					}
					if (area(i, j - 1) == areatype::INLET) {
						(*this)(i, j, 2, 0) = (*this)(i, j - 1, 4, 0) - 2.0 / 3.0 * (*this).density(i, j) * (c[4] * u_inlet);
					}
					if (area(i + 1, j) == areatype::INLET) {
						(*this)(i, j, 3, 0) = (*this)(i + 1, j, 1, 0) - 2.0 / 3.0 * (*this).density(i, j) * (c[1] * u_inlet);
					}
					if (area(i, j + 1) == areatype::INLET) {
						(*this)(i, j, 4, 0) = (*this)(i, j + 1, 2, 0) - 2.0 / 3.0 * (*this).density(i, j) * (c[2] * u_inlet);
					}
					if (area(i - 1, j - 1) == areatype::INLET) {
						(*this)(i, j, 5, 0) = (*this)(i - 1, j - 1, 7, 0) - 2.0 / 3.0 * (*this).density(i, j) * (c[7] * u_inlet);
					}
					if (area(i + 1, j - 1) == areatype::INLET) {
						(*this)(i, j, 6, 0) = (*this)(i + 1, j - 1, 8, 0) - 2.0 / 3.0 * (*this).density(i, j) * (c[8] * u_inlet);
					}
					if (area(i + 1, j + 1) == areatype::INLET) {
						(*this)(i, j, 7, 0) = (*this)(i + 1, j + 1, 5, 0) - 2.0 / 3.0 * (*this).density(i, j) * (c[5] * u_inlet);
					}
					if (area(i - 1, j + 1) == areatype::INLET) {
						(*this)(i, j, 8, 0) = (*this)(i - 1, j + 1, 6, 0) - 2.0 / 3.0 * (*this).density(i, j) * (c[6] * u_inlet);
					}

				}
			}
		}

		return *this;
	}

	template <int X, int Y>
	distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::oulet_outlflow_Neumann() {

		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (area(i, j) == areatype::outlet || area(i, j) == areatype::FB_outlet) {

					if (area(i - 1, j) == areatype::OUTLET) {
						(*this)(i, j, 1) = (*this)(i + 1, j, 1);
					}
					if (area(i, j - 1) == areatype::OUTLET) {
						(*this)(i, j, 2) = (*this)(i, j + 1, 2);
					}
					if (area(i + 1, j) == areatype::OUTLET) {
						(*this)(i, j, 3) = (*this)(i - 1, j, 3);
					}
					if (area(i, j + 1) == areatype::OUTLET) {
						(*this)(i, j, 4) = (*this)(i, j - 1, 4);
					}
					if (area(i - 1, j - 1) == areatype::OUTLET) {
						(*this)(i, j, 5) = (*this)(i + 1, j + 1, 5);
					}
					if (area(i + 1, j - 1) == areatype::OUTLET) {
						(*this)(i, j, 6) = (*this)(i - 1, j + 1, 6);
					}
					if (area(i + 1, j + 1) == areatype::OUTLET) {
						(*this)(i, j, 7) = (*this)(i - 1, j - 1, 7);
					}
					if (area(i - 1, j + 1) == areatype::OUTLET) {
						(*this)(i, j, 8) = (*this)(i + 1, j - 1, 8);
					}

				}
			}
		}

		return *this;
	}

	template <int X, int Y>
	distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::oulet_outlflow_convective() {

		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (area(i, j) == areatype::outlet || area(i, j) == areatype::FB_outlet) {
					//==========================================================================


				}
			}
		}

		return *this;
	}

}

//read area_file and design area
namespace distribution_function_template_space {

	//read area : S INLET OUTLET
	template<int X,int Y>
	int read_area_file(area_field<areatype,X,Y>& areafield,std::string& area_file_name) {

		std::fstream file(area_file_name, std::ios::in);

		if (!file) {
			std::cout << "open file error!" << std::endl;
		}

		std::vector<std::string> vec_str;
		std::vector<int> index;
		std::string str_for_cin;
		int count = 0;

		while (file >> str_for_cin) {
			vec_str.push_back(str_for_cin);
			if (str_for_cin == "block") {
				index.push_back(count);//0
			}
			else if (str_for_cin == "S") {
				index.push_back(count);//1
			}
			else if (str_for_cin == "INLET") {
				index.push_back(count);//2
			}
			else if (str_for_cin == "OUTLET") {
				index.push_back(count);//3
			}
			count++;
		}

		//detect
		if (index.size() != 4) {
			std::cout << "read area error!" << std::endl;
			return -1;
		}

		//int X,Y;
		//X = stoi(vec_str[index[0] + 1]);
		//Y = stoi(vec_str[index[0] + 2]);
		//std::cout << "X=" << X << " Y=" << Y << std::endl;

		//detect block
		//area_file 文件中"block"所包含的值与代码中的X Y对应，由于实际区域的长度为X+2 Y+2，所以S INLET OUTLET中的坐标上限为X+1或Y+1，下限为0
		if (stoi(vec_str[index[0] + 1]) != X && stoi(vec_str[index[0] + 2]) != Y ) {
			std::cout << "X or Y is not consistent with code" << std::endl;
		}


		//read S
		for (int i = index[1] + 1; i < index[2] - 1; i++, i++) {
			//std::cout << stoi(vec_str[i]) << " " << stoi(vec_str[i + 1]) << std::endl;
			areafield(stoi(vec_str[i]), stoi(vec_str[i + 1])) = areatype::S;
		}

		//read inlet
		for (int i = index[2] + 1; i < index[3] - 1; i++, i++) {
			//std::cout << stoi(vec_str[i]) << " " << stoi(vec_str[i + 1]) << std::endl;
			areafield(stoi(vec_str[i]), stoi(vec_str[i + 1])) = areatype::INLET;
		}

		//read outlet
		for (int i = index[3] + 1; i < count - 1; i++, i++) {
			//std::cout << stoi(vec_str[i]) << " " << stoi(vec_str[i + 1]) << std::endl;
			areafield(stoi(vec_str[i]), stoi(vec_str[i + 1])) = areatype::OUTLET;
		}
		
		file.close();
		return 0;
	}

	template<int X,int Y>
	void setF(area_field<areatype, X, Y>& areafield) {
		for (int i = 0; i <= X+1; i++) {
			for (int j = 0; j <= Y+1; j++) {
				if (areafield(i, j) != areatype::S && areafield(i, j) != areatype::INLET && areafield(i, j) != areatype::OUTLET) {
					areafield(i, j) = areatype::F;
				}
			}
		}
		return;
	}

	template<int X, int Y>
	void setFBFL(area_field<areatype, X, Y>& areafield) {
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (areafield(i, j) == areatype::F) {
					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (areafield(i + x, j + y) == areatype::S) {
								areafield(i, j) = areatype::FB;
							}
						}
					}
				}
			}
		}
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (areafield(i, j) == areatype::F) {
					areafield(i, j) = areatype::FL;
				}
			}
		}
		return;
	}

	template<int X, int Y>
	void setSBSL(area_field<areatype, X, Y>& areafield) {
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (areafield(i, j) == areatype::FB) {
					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (areafield(i + x, j + y) == areatype::S) {
								areafield(i + x, j + y) = areatype::SB;
							}
						}
					}
				}
			}
		}
		for (int i = 0; i <= X + 1; i++) {
			for (int j = 0; j <= Y + 1; j++) {
				if (areafield(i, j) == areatype::S) {
					areafield(i, j) = areatype::SL;
				}
			}
		}
		return;
	}

	template<int X, int Y>
	void setinlet_outlet(area_field<areatype, X, Y>& areafield) {
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (areafield(i, j) == areatype::FB || areafield(i, j) == areatype::FL){
					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (areafield(i + x, j + y) == areatype::INLET) {
								areafield(i + x, j + y) = areatype::inlet;
							}
							if (areafield(i + x, j + y) == areatype::OUTLET) {
								areafield(i + x, j + y) = areatype::outlet;
							}
						}
					}
				}
			}
		}
		return;
	}

	template<int X, int Y>
	void setFB_inlet_FB_outlet(area_field<areatype, X, Y>& areafield) {
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (areafield(i, j) == areatype::inlet) {
					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (areafield(i + x, j + y) == areatype::SB) {
								areafield(i + x, j + y) = areatype::FB_inlet;
							}
						}
					}
				}
				if (areafield(i, j) == areatype::outlet) {
					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (areafield(i + x, j + y) == areatype::SB) {
								areafield(i + x, j + y) = areatype::FB_outlet;
							}
						}
					}
				}
			}
		}
		return;
	}

	//inculded by tempalte class distribution_function_template_D2Q9<X, Y>
	template<int X, int Y>
	void distribution_function_template_D2Q9<X, Y>::setarea(std::string& area_file_name) {

		area_field<areatype, X, Y>& areafield = area;

		if(read_area_file<X, Y>(areafield, area_file_name) == -1) return;
		
		setF<X, Y>(areafield);

		setFBFL<X, Y>(areafield);

		//detect FB and FL
		for (int i = 0; i <= X + 1; i++) {
			for (int j = 0; j <= Y + 1; j++) {
				if (areafield(i, j) == areatype::FB || areafield(i, j) == areatype::FL) {
					if (i == 0 || i == X + 1 || j == 0 || j == Y + 1) {
						std::cout << "FB or FL is at margin e. i==0,j==0,i==X+1,j==Y+1 : error!" << std::endl;
						return;
					}
				}
			}
		}

		setSBSL<X, Y>(areafield);

		setinlet_outlet<X, Y>(areafield);

		//detect inlet and outlet
		for (int i = 0; i <= X + 1; i++) {
			for (int j = 0; j <= Y + 1; j++) {
				if (areafield(i, j) == areatype::inlet || areafield(i, j) == areatype::outlet) {
					if ((i == 1 || i == X) && (j == 1 || j == Y))
						;
					else {
						std::cout << "inlet or outlet is not near margin e. i!=1 and i!=X: error!" << std::endl;
						return;
					}
				}
			}
		}

		setFB_inlet_FB_outlet<X, Y>(areafield);

		return;
	}

}

//A distribution function class template for CGM.
namespace distribution_function_template_space {
	
	template <int X, int Y>
	class distribution_function_CGM_D2Q9 : public distribution_function_template_D2Q9<X,Y> {
	public:

		distribution_function_CGM_D2Q9(double rho_initial = 1,double nu_ = 1.0 / 6.0) :distribution_function_template_D2Q9<X,Y>(rho_initial, nu_){}

		using distribution_function_template_D2Q9<X, Y>::c;
		using distribution_function_template_D2Q9<X, Y>::w;
		using distribution_function_template_D2Q9<X, Y>::M;
		using distribution_function_template_D2Q9<X, Y>::InM;
		using distribution_function_template_D2Q9<X, Y>::S;
		using distribution_function_template_D2Q9<X, Y>::area;

		double& operator()(int i, int j, int q) {
			// (i,j,q) <-- [i-1][j-1][q] <-- (i-1) * Y * 9 + (j - 1) * 9 + q
			if (i > X) {
				printf("\ndistribution_function_CGM_D2Q9 is called, but the first subscript exceeds the limit: %d > X\n", i);
			}
			else if (i < 1) {
				printf("\ndistribution_function_CGM_D2Q9 is called, but the first subscript exceeds the limit: %d < 1\n", i);
			}
			else if (j > Y) {
				printf("\ndistribution_function_CGM_D2Q9 is called, but the second subscript exceeds the limit: %d > Y\n", j);
			}
			else if (j < 1) {
				printf("\ndistribution_function_CGM_D2Q9 is called, but the second subscript exceeds the limit: %d < 1\n", j);
			}
			else if (q < 0) {
				printf("\ndistribution_function_CGM_D2Q9 is called, but the third subscript exceeds the limit: %d < 0\n", q);
			}
			else if (q > 8) {
				printf("\ndistribution_function_CGM_D2Q9 is called, but the third subscript exceeds the limit: %d > 8\n", q);
			}
			int index = (i - 1) * Y * 9 + (j - 1) * 9 + q;
			return *(this->distribution_function_template_p + index);
		}

		//--1--SRT
		distribution_function_CGM_D2Q9<X, Y>& single_phase_collison_SRT(distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b);
		//--2--MRT
		distribution_function_CGM_D2Q9<X, Y>& single_phase_collison_MRT(distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b);

		//for color gradient model(CGM)
		static vector2D_field<X, Y> ns;
		static scalar_field_space::phase_field<X, Y> phaseField;
		static fluid_field<fluidtype, X, Y> fluidcolor;
		static vector<double, 9> B;
		static double contactAngle;
		static vector<double, 2> n_contact_angle;

		static double nu_is(int i, int j, distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b, double rho_r_reference = 1.0, double rho_b_reference = 1.0) {
			return (f_r.density(i, j) / rho_r_reference + f_b.density(i, j) / rho_b_reference) / (f_r.density(i, j) / rho_r_reference / f_r.nu + f_b.density(i, j) / rho_b_reference / f_b.nu);
		}

		static double omega_is(int i, int j, distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b, double rho_r_reference = 1.0, double rho_b_reference = 1.0) {
			return 2.0 / (6.0 * distribution_function_CGM_D2Q9<X, Y>::nu_is(i, j, f_r, f_b) + 1.0);
		}

		static void caculate_ns();

		distribution_function_CGM_D2Q9<X, Y>& wetting_boundary(double contact_angle, double delta_2 = 0.05);

		distribution_function_CGM_D2Q9<X, Y>& perturbationLiu2012(double A);

		//void perturbationLui2017();

		//This function is called by total distribution function.
		distribution_function_CGM_D2Q9<X, Y>& recolor(distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b,double beta);
		
		//This function is called by total distribution function.
		//This function is used in the color gradient model to get the total distribution function.
		distribution_function_CGM_D2Q9<X, Y>& blend(distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b);

	};

	//initialize
	//--1--
	template <int X, int Y>
	scalar_field_space::phase_field<X, Y> distribution_function_CGM_D2Q9<X, Y>::phaseField;
	//--2--
	template<int X, int Y>
	fluid_field<fluidtype, X, Y> distribution_function_CGM_D2Q9<X, Y>::fluidcolor(fluidtype::blue);
	//--3--
	template<int X, int Y>
	vector<double, 9> distribution_function_CGM_D2Q9<X, Y>::B({ -4.0 / 27.0,2.0 / 27.0,2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 5.0 / 108.0, 5.0 / 108.0 , 5.0 / 108.0 , 5.0 / 108.0 });
	//--4--
	template<int X, int Y>
	vector2D_field<X, Y> distribution_function_CGM_D2Q9<X, Y>::ns;
	//--5--
	template<int X, int Y>
	double distribution_function_CGM_D2Q9<X, Y>::contactAngle(90);
	//--6--
	template<int X, int Y>
	vector<double, 2> distribution_function_CGM_D2Q9<X, Y>::n_contact_angle({ cos(contactAngle / 180.0 * 3.141592653),sin(contactAngle / 180.0 * 3.141592653) });


	template <int X, int Y>
	inline distribution_function_CGM_D2Q9<X, Y>& distribution_function_CGM_D2Q9<X, Y>::blend(distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b) {
			for (int i = 1; i <= X; i++)
				for (int j = 1; j <= Y; j++)
					for (int q = 0; q <= 8; q++)
						(*this)(i, j, q) = f_r(i, j, q) + f_b(i, j, q);
		return *this;
	}

}

//1.collision 2.perturbation 3.recolor : the member functions OF distribution function tempalte class
namespace distribution_function_template_space {

	//SRT
	template <int X, int Y>
	distribution_function_CGM_D2Q9<X, Y>& distribution_function_CGM_D2Q9<X, Y>::single_phase_collison_SRT(distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b) {
		for (int i = 1; i <= X; i++)
			for (int j = 1; j <= Y; j++)
				for (int q = 0; q <= 8; q++)
					(*this)(i, j, q) = (1 - omega_is(i,j,f_r,f_b)) * (*this)(i, j, q) + omega_is(i, j, f_r, f_b) * (*this).equilibrium(i, j, q);
		return *this;
	}

	//MRT
	template <int X, int Y>
	distribution_function_CGM_D2Q9<X, Y>& distribution_function_CGM_D2Q9<X, Y>::single_phase_collison_MRT(distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b) {
		vector<double, 9> moment_pre_collison, moment_eq, moment_post_collison, f_vector_pre_collison, f_vector_eq, f_vector_post_collison;
		matrix<double, 9, 9> I(vector<double, 9>(1.0)), Matrix_S(S);
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				Matrix_S(8, 8) = Matrix_S(9, 9) = omega_is(i, j, f_r, f_b);
				for (int q = 0; q <= 8; q++)
					f_vector_pre_collison(q) = (*this)(i, j, q);

				for (int q = 0; q <= 8; q++)
					f_vector_eq(q) = (*this).equilibrium(i, j, q);

				//Convert to the moment space
				moment_pre_collison = M * f_vector_pre_collison;

				//Calculate the temporary equilibrium moment
				moment_eq = M * f_vector_eq;

				//Collision in moment space
				moment_post_collison = (I - Matrix_S) * moment_pre_collison + Matrix_S * moment_eq;

				//Back to particle space
				f_vector_post_collison = InM * moment_post_collison;

				for (int q = 0; q <= 8; q++)
					(*this)(i, j, q) = f_vector_post_collison(q);
			}
		}
		return *this;
	}

	template <int X, int Y>
	distribution_function_CGM_D2Q9<X, Y>& distribution_function_CGM_D2Q9<X, Y>::wetting_boundary(double contact_angle,double delta_2) {
		//wetting boundary conditions, which are used only for two-phase interfaces
		contactAngle = contact_angle;
		n_contact_angle = { cos(contact_angle / 180.0 * 3.141592653), sin(contact_angle / 180.0 * 3.141592653) };
		double temp1 = 0, temp2 = 0;
		//First step : estimate the phase field at SB
		for (int i = 0; i <= X + 1; i++) {
			for (int j = 0; j <= Y + 1; j++) {
				if (area(i,j) == areatype::SB) {
					temp1 = 0;
					temp2 = 0;
					phaseField(i, j) = 0;

					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (area(i + x, j + y) == areatype::FB) {
								temp1 += w(x * x + y * y) * phaseField(i + x, j + y, 0);
								temp2 += w(x * x + y * y);
							}
						}
					}

					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (area(i + x, j + y) == areatype::FB) {
								phaseField(i, j) = temp1 / temp2;
							}
						}
					}

				}

			}
		}
		//Second step : uses the phase field estimated at SB in the first step to caculate the phase field at FB
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (area(i, j) == areatype::FB) {

					phaseField.gradient(i, j)(0) = 1.0 / 12.0 * (phaseField(i + 1, j + 1) - phaseField(i - 1, j + 1))
						+ 1.0 / 3.0 * (phaseField(i + 1, j) - phaseField(i - 1, j)) + 1.0 / 12.0 * (phaseField(i + 1, j - 1) - phaseField(i - 1, j - 1));

					phaseField.gradient(i, j)(1) = 1.0 / 12.0 * (phaseField(i + 1, j + 1) - phaseField(i + 1, j - 1))
						+ 1.0 / 3.0 * (phaseField(i, j + 1) - phaseField(i, j - 1)) + 1.0 / 12.0 * (phaseField(i - 1, j + 1) - phaseField(i - 1, j - 1));

					if (phaseField.gradient(i, j) * phaseField.gradient(i, j) > delta_2) {
						fluidcolor = fluidtype::interface;
					}
					else if (phaseField(i,j) > 0) {
						fluidcolor = fluidtype::red;
					}
					else {
						fluidcolor = fluidtype::blue;
					}

				}
			}
		}
		//Third step : keep the magnitude of phase field gradient estimated at FB unchanged, then correct the direction of phase field gradient based the contact Angle and the fixed wall normal at FB
		double n1x, n1y, n2x, n2y, Fxn, Fyn, D1, D2, FF;
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (area(i,j) == areatype::FB) {

					if (fluidcolor(i,j) == fluidtype::interface) {

						n1x = vector<double, 2>({ ns(i, j)(0) ,-ns(i, j)(1) }) * n_contact_angle;
						n1y = vector<double, 2>({ ns(i, j)(1) ,ns(i, j)(0) }) * n_contact_angle;
						n2x = vector<double, 2>({ ns(i, j)(0) ,ns(i, j)(1) }) * n_contact_angle;
						n2y = vector<double, 2>({ ns(i, j)(1) ,-ns(i, j)(0) }) * n_contact_angle;

						FF = phaseField.gradient(i, j) * phaseField.gradient(i, j);
						Fxn = phaseField.gradient(i, j)(0) / sqrt(FF);
						Fyn = phaseField.gradient(i, j)(1) / sqrt(FF);
						D1 = sqrt((Fxn - n1x) * (Fxn - n1x) + (Fyn - n1y) * (Fyn - n1y));
						D2 = sqrt((Fxn - n2x) * (Fxn - n2x) + (Fyn - n2y) * (Fyn - n2y));
						if (D1 < D2) {
							phaseField.gradient(i, j)(0) = sqrt(FF) * n1x;
							phaseField.gradient(i, j)(1) = sqrt(FF) * n1y;
						}
						else if (D1 > D2) {
							phaseField.gradient(i, j)(0) = sqrt(FF) * n2x;
							phaseField.gradient(i, j)(1) = sqrt(FF) * n2y;
						}
						else {
							phaseField.gradient(i, j)(0) = sqrt(FF) * ns(i,j)(0);
							phaseField.gradient(i, j)(1) = sqrt(FF) * ns(i,j)(1);
						}
					}

				}
			}
		}
		return *this;
	}

	//perturbation
	template <int X, int Y>
	distribution_function_CGM_D2Q9<X, Y>& distribution_function_CGM_D2Q9<X, Y>::perturbationLiu2012(double A) {
		double FF = 0, Fc = 0;
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (fluidcolor(i,j) == fluidtype::interface) {

					FF = phaseField.gradient(i, j) * phaseField.gradient(i, j);

					for (int q = 0; q <= 8; q++) {
						Fc = phaseField.gradient(i, j) * c[q];
						(*this)(i, j, q) += A / 2.0 * sqrt(FF) * (w(q) * Fc * Fc / FF - B(q));
					}
				}
			}
		}
		return *this;
	}

	//recolor
	template <int X, int Y>
	distribution_function_CGM_D2Q9<X, Y>& distribution_function_CGM_D2Q9<X, Y>::recolor(distribution_function_CGM_D2Q9<X, Y>& f_r, distribution_function_CGM_D2Q9<X, Y>& f_b, double beta){
		double FF = 0, Fc = 0, rho, rho_r, rho_b;
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (fluidcolor(i, j) == fluidtype::interface) {

					FF = phaseField.gradient(i, j) * phaseField.gradient(i, j);
					rho = this->density(i, j);
					rho_r = f_r.density(i, j);
					rho_b = f_b.density(i, j);

					for (int q = 0; q <= 8; q++) {
						Fc = phaseField.gradient(i, j) * c[q];
						f_r(i, j, q) = rho_r / rho * (*this)(i, j, q) + beta * w(q) * rho_r * rho_b / rho *
							(c[q] * phaseField.gradient(i, j)) / sqrt((c[q] * c[q]) * (phaseField.gradient(i, j) * phaseField.gradient(i, j)));
						f_b(i, j, q) = rho_b / rho * (*this)(i, j, q) - beta * w(q) * rho_r * rho_b / rho *
							(c[q] * phaseField.gradient(i, j)) / sqrt((c[q] * c[q]) * (phaseField.gradient(i, j) * phaseField.gradient(i, j)));
					}
				}
			}
		}
		return *this;
	}

}

//design area the calculate ns
namespace distribution_function_template_space {
	
	

	double ww(int k) {
		if (k == 0)
			return 0;
		else if (k == 1)
			return 4.0 / 21.0;
		else if (k == 2)
			return 4.0 / 45.0;
		else if (k == 4)
			return 1.0 / 60.0;
		else if (k == 5)
			return 2.0 / 315.0;
		else if (k == 8)
			return 1.0 / 5040.0;
		return 0;
	}

	template<int X,int Y>
	double s(int& i, int& j, area_field<areatype, X, Y>& area) {
		
		if (i >= 0 && i <= X + 1 && j >= 0 && j <= Y + 1) {
			if (area(i, j) == areatype::SB || area(i, j) == areatype::SL)
				return 1.0;
			else
				return 0;
		}
		else
			return 1.0;

	}

	template<int X, int Y>
	double tempx(int& i, int& j, area_field<areatype, X, Y>& area) {
		double temp = 0;
		for (int x = -2; x <= 2; x++) {
			for (int y = -2; y <= 2; y++) {
				temp += ww(x * x + y * y) * s<x,y>(i + x, j + y, area) * x;
			}
		}
		return temp;
	}

	template<int X, int Y>
	double tempy(int& i, int& j, area_field<areatype, X, Y>& area) {
		double temp = 0;
		for (int x = -2; x <= 2; x++) {
			for (int y = -2; y <= 2; y++) {
				temp += ww(x * x + y * y) * s<X,Y>(i + x, j + y, area) * y;
			}
		}
		return temp;
	}

	template<int X,int Y>
	void distribution_function_CGM_D2Q9<X, Y>::caculate_ns() {

		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (area(i,j) == areatype::FB) {
					ns(i,j)(0) = tempx<X, Y>(i, j, area) / sqrt(tempx<X, Y>(i, j, area) * tempx<X, Y>(i, j, area) + tempy<X, Y>(i, j, area) * tempy<X, Y>(i, j, area));
					ns(i,j)(1) = tempy<X, Y>(i, j, area) / sqrt(tempx<X, Y>(i, j, area) * tempx<X, Y>(i, j, area) + tempy<X, Y>(i, j, area) * tempy<X, Y>(i, j, area));
				}
			}
		}
		//amend ns
		for (int i = 1; i <= X; i++) {
			for (int j = 1; j <= Y; j++) {
				if (area(i, j) != areatype::FB) {
					ns(i, j)(0) = 0;
					ns(i, j)(1) = 0;
				}
			}
		}

		return;
	}
}

#endif // !_DISTRIBUTION_FUNCTION_H_