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

		vector() :v(new T[N])
		{
			for (int i = 0; i < N; i++) {
				this->operator()(i) = 0;
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
				return *(vector2D_field<X, Y>::vector2D_field_p + index);
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
				return *(vector2D_field<X, Y>::vector2D_field_p + index);
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
			scalar_field(double scalar_initial = 0) : scalar_p(new double[X*Y])
			{
				for (int r = 0; r < X * Y; r++)
					*(scalar_p + r) = scalar_initial;
			}

			~scalar_field(){
				delete[]scalar_p;
			}

			double& operator()(int i,int j) {
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
			density_field(double density_initial = 1.0) :scalar_field<X, Y>(density_initial) {}

			double& operator()(int i, int j){
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

		};

		//define a phase_field template inheriting from scalar_field template
		template<int X, int Y>
		class phase_field :public scalar_field<X, Y> {
		public:

			phase_field(double phase_initial = 0) :scalar_field<X, Y>(phase_initial) {}

			double& operator()(int i, int j) {
				// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
				if (i > X) {
					printf("\nphase_field is called, but the first subscript exceeds the limit: %d > X\n", i);
				}
				else if (i < 1) {
					printf("\nphase_field is called, but the first subscript exceeds the limit: %d < 1\n", i);
				}
				else if (j > Y) {
					printf("\nphase_field is called, but the second subscript exceeds the limit: %d > Y\n", j);
				}
				else if (j < 1) {
					printf("\nphase_field is called, but the second subscript exceeds the limit: %d < 1\n", j);
				}
				int index = (i - 1) * Y + (j - 1);
				return *(scalar_field<X, Y>::scalar_p + index);
			}

			double& operator()(int i, int j,int) {
				// (i,j) <-- [i-1][j-1] <-- (i-1)*Y + (j-1)
				if (i > X) {
					printf("\nphase_field is called, but the first subscript exceeds the limit: %d > X\n", i);
				}
				else if (i < 1) {
					printf("\nphase_field is called, but the first subscript exceeds the limit: %d < 1\n", i);
				}
				else if (j > Y) {
					printf("\nphase_field is called, but the second subscript exceeds the limit: %d > Y\n", j);
				}
				else if (j < 1) {
					printf("\nphase_field is called, but the second subscript exceeds the limit: %d < 1\n", j);
				}
				int index = (i - 1) * Y + (j - 1);
				return *(scalar_field<X, Y>::scalar_p + index);
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
							printf("phase_field<X, Y>::calculate(rho1,rho2,1,1) is called,but at (%d,%d) rho1 + rho2 = %.5f< 0", i, j, b);
						}
					}
				}
				return *this;
			}

			//--2--phase_field's own function to calculate phase field and limit the scope to a certain fluid region F,FB,FL,inlet,outlet,inlet_S,outlet_S except S,SB,SL 
			//(r1-r2)/(r1+r2) or (r1/rr1-r2/rr2)/(r1/rr1+r2/rr2)
			phase_field<X, Y>& calculate(scalar_field_space::density_field<X, Y>& rho_1, scalar_field_space::density_field<X, Y>& rho_2, area_field<areatype,X,Y> areafield,double rho_reference_1 = 1, double rho_refernece_2 = 1) {
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
								printf("phase_field<X, Y>::calculate(rho1,rho2,areafield,1,1) is called,but at (%d,%d) rho1 + rho2 = %.5f< 0", i, j, b);
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
			phase_field<X, Y>& determine_fluidtype(fluid_field<fluidtype,X,Y> fluidfield, area_field<areatype,X,Y> areafield,double delta = 0.7){
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
								printf("phase_field<X, Y>::determine_fluidtype() is called,but ths value of phase_field at (%d,%d) is infinite.", i, j);
								return *this;
							}
						}
					}
				}
				return *this;
			}

			//--2--distinguish the fluidfield by the phase field gradient and delta
			//phase_field<X, Y>& determine_fluidtype(fluid_field<fluidtype, X, Y> fluidfield, area_field<areatype,X,Y> areafield,double delta_2 = 0.05) {
			//	
			//	for (int i = 1; i <= X; i++) {
			//		for (int j = 1; j <= Y; j++) {
			//			if (isfinite((*this)(i, j)) != 0) {

			//			}
			//		}
			//	}

			//	return *this;
			//}
			
			//solve for the phase field gradient
			vector2D_field<X, Y>& gradient(area_field<areatype, X, Y> areafield) {
				for (int i = 1; i <= X; i++) {
					for (int j = 1; j <= Y; j++) {
						if (areafield(i, j) != areatype::SL && areafield(i, j) != areatype::SB) {

							phase_field_gradient(i, j)(0) = 1.0 / 12.0 * ((*this)(i + 1, j + 1) - (*this)(i - 1, j + 1))
								+ 1.0 / 3.0 * ((*this)(i + 1, j) - (*this)(i - 1, j)) + 1.0 / 12.0 * ((*this)(i + 1, j - 1) - (*this)(i - 1, j - 1));

							phase_field_gradient(i, j)(1) = 1.0 / 12.0 * ((*this)(i + 1, j + 1) - (*this)(i + 1, j - 1))
								+ 1.0 / 3.0 * ((*this)(i, j + 1) - (*this)(i, j - 1)) + 1.0 / 12.0 * ((*this)(i - 1, j + 1) - (*this)(i - 1, j - 1));
							
						}
					}
				}
				return phase_field_gradient;
			}

			protected:
				//phase field gradient
				vector2D_field<X, Y> phase_field_gradient;

		};

	}

	class phase_field_gradient {
	public:
		phase_field_gradient() = default;
		phase_field_gradient(int n) :phase_field_gradient_site(new vector<double,2>[n]){

		}
		~phase_field_gradient()
		{
			delete[]phase_field_gradient_site;
		}
		vector<double, 2> *phase_field_gradient_site;
	};

}

//define classes associated with area type and fluid type that represents different areas
namespace distribution_function_template_space {

	enum class areatype{S,F,SB,SL,FB,FL,inlet,outlet,inlet_S,outlet_S};
	enum class fluidtype{red,blue,interface};

	template<typename T,int X, int Y>
	class area_field {
	public:
		area_field(T area_initial_type) :area_field_p(new T[(X + 2) * (Y + 2)]) {
			for (int r = 0; r < (X + 2) * (Y + 2); r++) {
				*area_field_p = area_initial_type;
			}
		}
		~area_field() {
			delete[]area_field_p;
		}
		T& operator()(int i,int j) {
			// (i,j) <-- [i][j] <-- i * Y + j
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
			int index = i * Y + j;
			return *(area_field_p + index);
		}
	protected:
		T* area_field_p;
	};

	template<typename T, int X, int Y>
	class fluid_field {
	public:
		fluid_field(T fluid_initial_type) :fluid_field_p(new T[X * Y]) {
			for (int r = 0; r < X * Y; r++) {
				*fluid_field_p = fluid_initial_type;
			}
		}
		~fluid_field() {
			delete[]fluid_field_p;
		}
		T& operator()(int i, int j) {
			// (i,j) <-- [i-1][j-1] <-- (i-1) * Y + (j-1)
			if (i > X) {
				printf("\nfluid_field is called, but the first subscript exceeds the limit: %d > X + 1\n", i);
			}
			else if (i < 1) {
				printf("\nfluid_field is called, but the first subscript exceeds the limit: %d < 0\n", i);
			}
			else if (j > Y) {
				printf("\nfluid_field is called, but the second subscript exceeds the limit: %d > Y + 1\n", j);
			}
			else if (j < 1) {
				printf("\nfluid_field is called, but the second subscript exceeds the limit: %d < 0\n", j);
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
		
		distribution_function_template_D2Q9<X, Y>& blend(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2); //This function is used in the color gradient model to get the total distribution function.
		
		distribution_function_template_D2Q9<X, Y>& streaming();

		//Solve for the equilibrium distribution function.
		//distribution_function_template_D2Q9.equilibrium(velocity2D_field,scalar_field);
		distribution_function_template_D2Q9<X,Y>& equilibrium(vector2D_field_space::velocity2D_field<X,Y> &velocity,scalar_field_space::scalar_field<X,Y> &rho);

		//Detect if the value of the distribution function is abnormal.
		bool detect();

	protected:

		distribution_function_template_space::vector<double, 9> w;
		distribution_function_template_space::vector<double, 2> c[9];
		double* distribution_function_template_p;
		
	};

	template <int X, int Y>
	inline distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::streaming() {
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
		return *this;
	}

	//template <int X,int Y>
	//distribution_function_template_D2Q9<X, Y> operator+(distribution_function_template_D2Q9<X, Y>& f1, distribution_function_template_D2Q9<X, Y>& f2) {
	//}

	template <int X,int Y>
	inline distribution_function_template_D2Q9<X, Y>& distribution_function_template_D2Q9<X, Y>::blend(distribution_function_template_D2Q9<X,Y>& f1, distribution_function_template_D2Q9<X, Y>& f2) {
		for (int q = 0; q <= 8; q++)
			for (int i = 1; i <= X; i++)
				for (int j = 1; j <= Y; j++)
					(*this)(i, j, q) = f1(i, j, q) + f2(i, j, q);
		return *this;
	}

	template <int X, int Y>
	inline distribution_function_template_D2Q9<X,Y>& distribution_function_template_D2Q9<X, Y>::equilibrium( vector2D_field_space::velocity2D_field<X, Y>& velocity,
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

//1.collision 2.perturbation 3.recolor
namespace distribution_function_template_space {

}

#endif // !_DISTRIBUTION_FUNCTION_H_
