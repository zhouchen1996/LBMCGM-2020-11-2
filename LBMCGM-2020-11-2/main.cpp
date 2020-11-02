#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "MRT.h"
//��������"Lattice Boltzmann simulation of pressure-driven two-phase flows in capillary tube and porous medium"
//�Լ�huang2014�����¶Գ���ڱ߽���е���
//�������������(��ɫ):���Ϊ�ٶȱ߽磬����Ϊѹ��Ϊ0�ı߽�
//���ڿ����������(��ɫ):���Ϊѹ��Ϊ0�ı߽磬����Ϊ��ѹ�ı߽�
//������ԣ�   
//		��ɫ :	��� : �ٶ�Ϊ u0
//				���� : ѹ��Ϊ ��Сֵ ����
//		��ɫ :  �������ھ�Ϊѹ����Сֵ
//�޸���sourceTermPrecondition()�����������ڽ���ı�����������ֵΪ0

//��������"Lattice Boltzmann simulation of pressure-driven two-phase flows in capillary tube and porous medium"
//���ֲܷ������ĽǶ������ڳ��ڱ߽�

using namespace std;
constexpr int nx = 100;
constexpr int ny = 30;
constexpr int Q = 9;
double N_r[nx][ny][Q]{ 0 };
double N_b[nx][ny][Q]{ 0 };
double N[nx][ny][Q]{ 0 };
double Neq_r[nx][ny][Q]{ 0 };
double Neq_b[nx][ny][Q]{ 0 };
double rho_r[nx][ny]{ 0 };
double rho_b[nx][ny]{ 0 };
double rho[nx][ny]{ 0 };
double rhoIn = 1e-3;//����Ϊ��ɫ����׼�������0ѹ������ֵ
double rhoOut = 0.1;//�����ڳ��ڴ���������ܶ�
double u[nx][ny]{ 0 };
double v[nx][ny]{ 0 };
double uInlet = 1e-3;//�������
double uInletDistribution[ny];
double vInlet = 0.0;
//ע�⣬����uInlet rhoOutl�Ĵ�С�Ա�֤���������0.1���¡���ǰ1000��������0.1���²��ܱ�֤�������ȶ���
double Phi[nx][ny]{ 0 };//�ೡ
double F_x[nx][ny]{ 0 };//�ೡ�ݶ�
double F_y[nx][ny]{ 0 };
double n[nx][ny][2]{ 0 };//�ೡ�ݶȵĵ�λ����
double kappa[nx][ny]{ 0 };//����
double F[nx][ny][2]{ 0 };//������������������������������Ҳ���Լ�����������
double cx[9]{ 0,1,0,-1,0,1,-1,-1,1 };
double cy[9]{ 0,0,1,0,-1,1,1,-1,-1 };
double w[9] = { 4.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0 };
double rho_r0 = 1;//��ʼ�ܶ�
double rho_b0 = 1;
double nu_r = 1.0 / 6.0;//ճ��ϵ����ճ��Խ�󣬳�ԥ�ʾ�Խ��,����ƽ����ٶ�Ҳ��Խ��,���Ҳ��Ӱ��α���������������ճ�Զ����ʱ��α�������С
double nu_b = 1.0 / 60.0;
double contactAngle = 120.0 / 180.0 * 3.141592653;//�Ӵ���!!!!!!!!!!!!!!!!!
double theta = 0.004;//��������,�����������ϴ�ʱ���߽�᲻�ȶ�,���Ǹ�����Ҫ�Ĳ��������������������ƣ�α����Ҳ����Ӱ�죬�������ԣ���ò�����0.005
double beta = 0.99;
double delta2 = 0.05;//�����ೡ�ݶȿ��ƽ�����,Ŀǰ���ֵ��ȡ���Ǹ���ͼ���Գ�����
double uu, cu, FF, Fc;
// S�ԽǾ��� �ɳ�ϵ������ s7��s8��BGK�е�omega��ͬ 
double s0 = 1.0, s1 = 1.63, s2 = 1.14, s3 = 1.0, s4 = 1.92, s5 = 1.0, s6 = 1.92, s7, s8;//s7 s8��仯
double S[9] = { s0,s1,s2,s3,s4,s5,s6,s7,s8 };
// ƽ��طֲ�����
double temp_moment[9]{ 0 };
double temp_moment_eq[9]{ 0 };
double temp_moment_post_collison[9]{ 0 };
double temp_moment_F[9]{ 0 };
double temp_moment_post_collison_F[9]{ 0 };
enum areatype { interface, red, blue, other };
areatype area[nx][ny];

int interval = 10;

namespace Region {
	enum regiontype { F, FB, FL, S, SB, SL };//һ����FB,FL,SB,SL
	regiontype region[nx][ny];
	double ns[nx][ny][2];
	double w(int k) {
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
	double s(int i, int j) {
		if (region[i][j] == SB || region[i][j] == SL) {
			return 1.0;
		}
		else if (region[i][j] == FB || region[i][j] == FL) {
			return 0.0;
		}
		else return 0;
	}
	double tempx(int& i, int& j) {
		double temp = 0;
		for (int x = -2; x <= 2; x++) {
			for (int y = -2; y <= 2; y++) {
				temp += w(x * x + y * y) * s(i + x, j + y) * x;
			}
		}
		return temp;
	}
	double tempy(int& i, int& j) {
		double temp = 0;
		for (int x = -2; x <= 2; x++) {
			for (int y = -2; y <= 2; y++) {
				temp += w(x * x + y * y) * s(i + x, j + y) * y;
			}
		}
		return temp;
	}
	void cal_ns() {
		for (int i = 2; i < nx - 2; i++) {
			for (int j = 2; j < ny - 2; j++) {
				if (region[i][j] == FB) {
					ns[i][j][0] = tempx(i, j) / sqrt(tempx(i, j) * tempx(i, j) + tempy(i, j) * tempy(i, j));
					ns[i][j][1] = tempy(i, j) / sqrt(tempx(i, j) * tempx(i, j) + tempy(i, j) * tempy(i, j));
				}
			}
		}
		//����ns
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (region[i][j] != FB) {
					ns[i][j][0] = 0;
					ns[i][j][1] = 0;
				}
			}
		}
		return;
	}
}

namespace Region {
	namespace setRegion {
		//��setS()������ƺù��������ֱ�ӵ���Region::setRegion::setProcess();
		void setS() {
			//ֻҪ��ƺ���������Ϳ�����

			double length;
			double obstacle_cylinder_centerx[2] = { nx / 3.0, 2.0 * nx / 3.0 };
			double obstacle_cylinder_certery[2] = { ny / 2.0, ny / 2.0 };
			double obstacle_cylinder_r = ny / 4.0;

			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					if (j == 0 || j == 1 || j == ny - 1 || j == ny - 2) {
						region[i][j] = S;
					}

					for (int k = 0; k < 2; k++) {
						length = sqrt(double(pow(i - obstacle_cylinder_centerx[k], 2) + pow(j - obstacle_cylinder_certery[k], 2)));
						if (length < obstacle_cylinder_r)
						{
							region[i][j] = S;
						}
					}

				}
			}
			return;
		}
		void setF() {
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					if (region[i][j] != S)
						region[i][j] = F;
				}
			}
			return;
		}
		void setFBFL() {
			for (int i = 1; i < nx - 1; i++) {
				for (int j = 1; j < ny - 1; j++) {
					if (region[i][j] == F) {
						for (int x = -1; x <= 1; x++) {
							for (int y = -1; y <= 1; y++) {
								if (region[i + x][j + y] == S) {
									region[i][j] = FB;
								}
							}
						}
					}
				}
			}
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					if (region[i][j] == F) {
						region[i][j] = FL;
					}
				}
			}
			return;
		}
		void setSBSL() {
			for (int i = 1; i < nx - 1; i++) {
				for (int j = 1; j < ny - 1; j++) {
					if (region[i][j] == S) {
						for (int x = -1; x <= 1; x++) {
							for (int y = -1; y <= 1; y++) {
								if (region[i + x][j + y] == FB) {
									region[i][j] = SB;
								}
							}
						}
					}
				}
			}
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					if (region[i][j] == S) {
						region[i][j] = SL;
					}
				}
			}
			return;
		}
		void setProcess() {
			setS();
			setF();
			setFBFL();
			setSBSL();
			return;
		}
	}
}

namespace wettingBoundary {

	double w(int k) {
		if (k == 0)
			return 4.0 / 9.0;
		else if (k == 1)
			return 1.0 / 9.0;
		else if (k == 2)
			return 1.0 / 36.0;
		else return 0;
	}
	void wetting_boundary(double(&F_x)[nx][ny], double(&F_y)[nx][ny]) {
		//ʪ��߽�����,����߽�����ֻʹ�����������
		//�ڱ߽�ʹ��������ݶȣ�ֻ��ͬʱ�ڽӴ��ߺͱ߽��ϵľ�������
		double temp1, temp2;
		//��һ���ȹ���SB����S���ೡPhi
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				if (Region::region[i][j] == Region::SB) {
					temp1 = 0;
					temp2 = 0;
					Phi[i][j] = 0;

					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (Region::region[i + x][j + y] == Region::FB) {
								temp1 += w(x * x + y * y) * Phi[i + x][j + y];
								temp2 += w(x * x + y * y);
							}
						}
					}

					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (Region::region[i + x][j + y] == Region::FB) {
								Phi[i][j] = temp1 / temp2;
							}
						}
					}

				}

			}
		}
		//�ڶ������õ�һ�������SB�����ೡPhi����FB�����ೡ�ݶ�F_x,F_y
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				if (Region::region[i][j] == Region::FB) {

					F_x[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i - 1][j + 1])
						+ 2.0 / 3.0 * (Phi[i + 1][j] - Phi[i - 1][j]) + 1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i - 1][j - 1]));

					F_y[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i + 1][j - 1])
						+ 2.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i - 1][j + 1] - Phi[i - 1][j - 1]));

					//��������2 ���ೡ�ݶȵĴ�С����aqrt(delta2)ʱ��Ϊ����
					if (F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j] > delta2) {
						area[i][j] = interface;
					}
					else {
						area[i][j] = other;
					}
				}
			}
		}
		//������������FB������ೡ�ݶ�F_x,F_y��С���䣬���ݽӴ�����FB���Ḻ́ڷ���ns��������
		double n1x, n1y, n2x, n2y, Fxn, Fyn, D1, D2;
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				if (Region::region[i][j] == Region::FB) {

					if (area[i][j] == interface) {
						n1x = Region::ns[i][j][0] * cos(contactAngle) - Region::ns[i][j][1] * sin(contactAngle);
						n1y = Region::ns[i][j][1] * cos(contactAngle) + Region::ns[i][j][0] * sin(contactAngle);
						n2x = Region::ns[i][j][0] * cos(contactAngle) + Region::ns[i][j][1] * sin(contactAngle);
						n2y = Region::ns[i][j][1] * cos(contactAngle) - Region::ns[i][j][0] * sin(contactAngle);
						FF = F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j];
						Fxn = F_x[i][j] / sqrt(FF);
						Fyn = F_y[i][j] / sqrt(FF);
						D1 = sqrt((Fxn - n1x) * (Fxn - n1x) + (Fyn - n1y) * (Fyn - n1y));
						D2 = sqrt((Fxn - n2x) * (Fxn - n2x) + (Fyn - n2y) * (Fyn - n2y));
						if (D1 < D2) {
							F_x[i][j] = sqrt(FF) * n1x;
							F_y[i][j] = sqrt(FF) * n1y;
						}
						else if (D1 > D2) {
							F_x[i][j] = sqrt(FF) * n2x;
							F_y[i][j] = sqrt(FF) * n2y;
						}
						else {
							F_x[i][j] = sqrt(FF) * Region::ns[i][j][0];
							F_y[i][j] = sqrt(FF) * Region::ns[i][j][1];
						}
					}

				}
			}
		}
		return;
	}
}

namespace inoutBoundary {

	void caculate_uInletDistribution() {
		uInletDistribution[0] = 0;
		uInletDistribution[1] = 0;
		uInletDistribution[ny - 1] = 0;
		uInletDistribution[ny - 2] = 0;
		for (int j = 2; j <= ny - 3; j++) {
			uInletDistribution[j] = -4.0 * uInlet / ((ny - 5.0) * (ny - 5.0)) * (j - 2.0) * (3.0 - ny + j);
		}
		return;
	}

	void inletVelocity(double(&f)[nx][ny][Q], double(&f_r)[nx][ny][Q], double(&f_b)[nx][ny][Q],
		double(&rho)[nx][ny], double(&rho_r)[nx][ny], double(&rho_b)[nx][ny],
		double(&ux)[nx][ny], double(&uy)[nx][ny], double(&u)[ny]) {

		//�ϲ��ֲ�����
		for (int j = 2; j <= ny - 3; j++) {
			for (int k = 0; k < Q; k++) {
				f[0][j][k] = f_r[0][j][k] + f_b[0][j][k];
			}
		}

		//��ڵĵ�λ����Ϊ(1,0)ʱ���ٶȱ߽磬(1,0)ָ������
		for (int j = 3; j <= ny - 4; j++) {
			rho[0][j] = ((f[0][j][0] + f[0][j][2] + f[0][j][4]) + 2.0 * (f[0][j][3] + f[0][j][6] + f[0][j][7])) / (1 - u[j]);
			f[0][j][1] = f[0][j][3] + 2.0 * rho[0][j] * u[j] / 3.0;
			f[0][j][5] = f[0][j][7] - 0.5 * (f[0][j][2] - f[0][j][4]) + rho[0][j] * u[j] / 6.0 + 0.5 * rho[0][j] * vInlet;
			f[0][j][8] = f[0][j][6] + 0.5 * (f[0][j][2] - f[0][j][4]) + rho[0][j] * u[j] / 6.0 - 0.5 * rho[0][j] * vInlet;
		}

		//��ڵ�corner����Ҫ�ر���
		//���½�
		f[0][2][1] = f[0][2][3]; f[0][2][2] = f[0][2][4]; f[0][2][5] = f[0][2][7];
		f[0][2][6] = f[0][2][8] = 0.5 * (rho[0][3] - (f[0][2][0] + f[0][2][1] + f[0][2][2] + f[0][2][3] + f[0][2][4] + f[0][2][5] + f[0][2][7]));
		//���Ͻ�
		f[0][ny - 3][1] = f[0][ny - 3][3]; f[0][ny - 3][4] = f[0][ny - 3][2]; f[0][ny - 3][8] = f[0][ny - 3][6];
		f[0][ny - 3][5] = f[0][ny - 3][7] = 0.5 * (rho[0][ny - 4] - (f[0][ny - 3][0] + f[0][ny - 3][1] + f[0][ny - 3][2] + f[0][ny - 3][3] + f[0][ny - 3][4] + f[0][ny - 3][6] + f[0][ny - 3][8]));

		//�����ܶȷֲ��ֽ�ֲ�����
		for (int j = 2; j <= ny - 3; j++) {
			for (int k = 0; k < Q; k++) {
				f_r[0][j][k] = f[0][j][k] * rho_r[0][j] / rho[0][j];
				f_b[0][j][k] = f[0][j][k] * rho_b[0][j] / rho[0][j];
			}
		}

		return;
	}

	void inletPressure(double(&f)[nx][ny][Q], double(&f_r)[nx][ny][Q], double(&f_b)[nx][ny][Q],
		double(&rho)[nx][ny], double(&rho_r)[nx][ny], double(&rho_b)[nx][ny],
		double(&ux)[nx][ny], double(&uy)[nx][ny], double& rhoIn) {
		//��������ֱ���ò�ͬ������ܶ�
		//��ڵĵ�λ����Ϊ(1,0)ʱ��ѹ���߽�,(1,0)ָ������

		//�ϲ��ֲ�����
		for (int j = 2; j <= ny - 3; j++) {
			for (int k = 0; k < Q; k++) {
				f[0][j][k] = f_r[0][j][k] + f_b[0][j][k];
			}
		}

		for (int j = 3; j <= ny - 4; j++) {
			u[0][j] = 1.0 - 1.0 / rhoIn * (2 * (f[0][j][3] + f[0][j][6] + f[0][j][7]) +
				(f[0][j][0] + f[0][j][2] + f[0][j][4]));
			f[0][j][1] = f[0][j][3] + 2.0 / 3.0 * rhoIn * u[0][j];
			f[0][j][5] = f[0][j][7] - 0.5 * (f[0][j][2] - f[0][j][4]) + 1.0 / 6.0 * rhoIn * u[0][j];
			f[0][j][8] = f[0][j][6] + 0.5 * (f[0][j][2] - f[0][j][4]) + 1.0 / 6.0 * rhoIn * u[0][j];
		}

		//���ڵ�corner����Ҫ�ر���
		//���Ͻ�
		f[0][ny - 3][1] = f[0][ny - 3][3];
		f[0][ny - 3][4] = f[0][ny - 3][2];
		f[0][ny - 3][8] = f[0][ny - 3][6];
		f[0][ny - 3][5] = f[0][ny - 3][7] = 0.5 * (rhoIn - (f[0][ny - 3][1] + f[0][ny - 3][2] +
			f[0][ny - 3][3] + f[0][ny - 3][4] + f[0][ny - 3][6] + f[0][ny - 3][8]));
		//���½�
		f[0][2][1] = f[0][2][3];
		f[0][2][5] = f[0][2][7];
		f[0][2][4] = f[0][2][2];
		f[0][2][6] = f[0][2][8] = 0.5 * (rhoIn - (f[0][2][1] + f[0][2][2] +
			f[0][2][3] + f[0][2][4] + f[0][2][5] + f[0][2][7]));

		//�����ܶȷֲ�����������ķֲ�����
		for (int j = 2; j <= ny - 3; j++) {
			for (int k = 0; k < Q; k++) {
				f_r[0][j][k] = f[0][j][k] * rho_r[0][j] / rho[0][j];
				f_b[0][j][k] = f[0][j][k] * rho_b[0][j] / rho[0][j];
			}
		}

		return;
	}

	void outletPressure(double(&f)[nx][ny][Q], double(&f_r)[nx][ny][Q], double(&f_b)[nx][ny][Q],
		double(&rho)[nx][ny], double(&rho_r)[nx][ny], double(&rho_b)[nx][ny],
		double(&ux)[nx][ny], double(&uy)[nx][ny], double& rhoOut) {
		//��������ֱ���ò�ͬ�ĳ����ܶ�
		//���ڵĵ�λ����Ϊ(1,0)ʱ��ѹ���߽�,(1,0)ָ������

		//�ϲ��ֲ�����
		for (int j = 2; j <= ny - 3; j++) {
			for (int k = 0; k < Q; k++) {
				f[nx - 1][j][k] = f_r[nx - 1][j][k] + f_b[nx - 1][j][k];
			}
		}

		for (int j = 3; j <= ny - 4; j++) {
			u[nx - 1][j] = 1.0 / rhoOut * (2 * (f[nx - 1][j][1] + f[nx - 1][j][5] + f[nx - 1][j][8]) +
				(f[nx - 1][j][0] + f[nx - 1][j][2] + f[nx - 1][j][4])) - 1.0;
			f[nx - 1][j][3] = f[nx - 1][j][1] - 2.0 / 3.0 * rhoOut * u[nx - 1][j];
			f[nx - 1][j][6] = f[nx - 1][j][8] - 0.5 * (f[nx - 1][j][2] - f[nx - 1][j][4]) - 1.0 / 6.0 * rhoOut * u[nx - 1][j];
			f[nx - 1][j][7] = f[nx - 1][j][5] + 0.5 * (f[nx - 1][j][2] - f[nx - 1][j][4]) - 1.0 / 6.0 * rhoOut * u[nx - 1][j];
		}

		//���ڵ�corner����Ҫ�ر���
		//���Ͻ�
		f[nx - 1][ny - 3][3] = f[nx - 1][ny - 3][1];
		f[nx - 1][ny - 3][7] = f[nx - 1][ny - 3][5];
		f[nx - 1][ny - 3][4] = f[nx - 1][ny - 3][2];
		f[nx - 1][ny - 3][6] = f[nx - 1][ny - 3][8] = 0.5 * (rhoOut - (f[nx - 1][ny - 3][1] + f[nx - 1][ny - 3][2] +
			f[nx - 1][ny - 3][3] + f[nx - 1][ny - 3][4] + f[nx - 1][ny - 3][5] + f[nx - 1][ny - 3][7]));
		//���½�
		f[nx - 1][2][2] = f[nx - 1][2][4];
		f[nx - 1][2][3] = f[nx - 1][2][1];
		f[nx - 1][2][6] = f[nx - 1][2][8];
		f[nx - 1][2][5] = f[nx - 1][2][7] = 0.5 * (rhoOut - (f[nx - 1][2][1] + f[nx - 1][2][2] +
			f[nx - 1][2][3] + f[nx - 1][2][4] + f[nx - 1][2][6] + f[nx - 1][2][8]));

		//�����ܶȷֲ�����������ķֲ�����
		for (int j = 2; j <= ny - 3; j++) {
			for (int k = 0; k < Q; k++) {
				f_r[nx - 1][j][k] = f[nx - 1][j][k] * rho_r[nx - 1][j] / rho[nx - 1][j];;
				f_b[nx - 1][j][k] = f[nx - 1][j][k] * rho_b[nx - 1][j] / rho[nx - 1][j];;
			}
		}

		return;
	}

}

double nu_(int& i, int j) {
	return (rho_r[i][j] + rho_b[i][j]) / (rho_r[i][j] / nu_r + rho_b[i][j] / nu_b);
}

double omega(int& i, int& j) {
	return 2.0 / (6.0 * nu_(i, j) + 1.0);
}

double cosphi(int& k, double& Fc, double& FF) {
	if (k == 0) {
		return 0;
	}
	else if (k == 1 || k == 2 || k == 3 || k == 4) {
		double temp1 = Fc / sqrt(FF);
		return ((isnormal(temp1) == 1) ? temp1 : 0);
	}
	else {
		double temp2 = Fc / sqrt(FF * 2.0);
		return ((isnormal(temp2) == 1) ? temp2 : 0);
	}
}

double cosphiImprove(int& k, double& Fc, double& FF) {
	if (k == 0) {
		return 0;
	}
	else if (k == 1 || k == 2 || k == 3 || k == 4) {
		double temp1 = Fc / sqrt(FF);
		return ((isnormal(temp1) == 1) ? temp1 : 0);
	}
	else {
		double temp2 = Fc / sqrt(FF);
		return ((isnormal(temp2) == 1) ? temp2 : 0);
	}
}

void streaming(double(&f)[nx][ny][Q]) {
	//������
	//����Ϊ streaming(N_r);
	for (int j = 0; j < ny; j++) {
		for (int i = nx - 1; i > 0; i--) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB || Region::region[i][j] == Region::SB) {
				f[i][j][1] = f[i - 1][j][1];
			}
		}
	}
	for (int i = 0; i < nx; i++) {
		for (int j = ny - 1; j > 0; j--) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB || Region::region[i][j] == Region::SB) {
				f[i][j][2] = f[i][j - 1][2];
			}
		}
	}
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx - 1; i++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB || Region::region[i][j] == Region::SB) {
				f[i][j][3] = f[i + 1][j][3];
			}
		}
	}
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny - 1; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB || Region::region[i][j] == Region::SB) {
				f[i][j][4] = f[i][j + 1][4];
			}
		}
	}
	for (int j = ny - 1; j > 0; j--) {
		for (int i = nx - 1; i > 0; i--) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB || Region::region[i][j] == Region::SB) {
				f[i][j][5] = f[i - 1][j - 1][5];
			}
		}
	}
	for (int j = ny - 1; j > 0; j--) {
		for (int i = 0; i < nx - 1; i++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB || Region::region[i][j] == Region::SB) {
				f[i][j][6] = f[i + 1][j - 1][6];
			}
		}
	}
	for (int j = 0; j < ny - 1; j++) {
		for (int i = 0; i < nx - 1; i++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB || Region::region[i][j] == Region::SB) {
				f[i][j][7] = f[i + 1][j + 1][7];
			}
		}
	}
	for (int j = 0; j < ny - 1; j++) {
		for (int i = nx - 1; i > 0; i--) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB || Region::region[i][j] == Region::SB) {
				f[i][j][8] = f[i - 1][j + 1][8];
			}
		}
	}
	return;
}

void out_file(const string& str,			//strΪ�ļ���
	const string& scalar_attribute_1,		//scalar_attribute_1Ϊ�������Ե���������
	const string& vector_attribute_1,		//vector_attribute_1Ϊʸ�����Ե���������
	const double(&u_vector_x)[nx][ny],		//u_vector_xΪ�ٶȳ���x����ķ���
	const double(&v_vector_y)[nx][ny],		//v_vector_yΪ�ٶȳ���y����ķ���
	const double(&Phi)[nx][ny]) {			//rho_scalar_1Ϊ�ೡ
	//�����ļ��������
	//���Ƕ������ĺ���������ļ���ʽ
	//����Ϊ��
	/*
	if (step % interval == 0) {
		out_file("Phi" + to_string(step / interval) + ".vtk", "Phi", "velocity", u, v, Phi);
	}
	*/
	ofstream outfile;
	outfile.open(str);
	outfile << "# vtk DataFile Version 3.0" << endl;
	outfile << str << endl;
	outfile << "ASCII" << endl;
	outfile << "DATASET STRUCTURED_GRID" << endl;
	outfile << "DIMENSIONS" << " " << nx << " " << ny << " " << 1 << endl;
	outfile << "POINTS" << " " << nx * ny * 1 << " " << "double" << endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			outfile << i << " " << j << " " << 0 << endl;
		}
	}
	outfile << "POINT_DATA" << " " << nx * ny * 1 << endl;
	outfile << "SCALARS" << " " << scalar_attribute_1 << " " << "double 1" << endl;
	outfile << "LOOKUP_TABLE default" << endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			outfile << Phi[i][j] << endl;
		}
	}
	outfile << "VECTORS" << " " << vector_attribute_1 << " " << "float" << endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			outfile << u_vector_x[i][j] << " " << v_vector_y[i][j] << " " << 0 << endl;
		}
	}
	outfile.close();
	return;
}

void out_file(const string& str,			//strΪ�ļ���
	const string& scalar_attribute_1,		//scalar_attribute_1Ϊ�������Ե���������
	const string& scalar_attribute_2,
	const string& vector_attribute_1,		//vector_attribute_1Ϊʸ�����Ե���������
	const double(&u_vector_x)[nx][ny],		//u_vector_xΪ�ٶȳ���x����ķ���
	const double(&v_vector_y)[nx][ny],		//v_vector_yΪ�ٶȳ���y����ķ���
	const double(&Phi)[nx][ny],
	const areatype(&area)[nx][ny]) {
	//�����ļ��������
	//���Ƕ������ĺ���������ļ���ʽ
	//����Ϊ��
	/*
	if (step % interval == 0) {
		out_file("Phi" + to_string(step / interval) + ".vtk", "Phi", "area","velocity", u, v, Phi,area);
	}
	*/
	ofstream outfile;
	outfile.open(str);
	outfile << "# vtk DataFile Version 3.0" << endl;
	outfile << str << endl;
	outfile << "ASCII" << endl;
	outfile << "DATASET STRUCTURED_GRID" << endl;
	outfile << "DIMENSIONS" << " " << nx << " " << ny << " " << 1 << endl;
	outfile << "POINTS" << " " << nx * ny * 1 << " " << "double" << endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			outfile << i << " " << j << " " << 0 << endl;
		}
	}
	outfile << "POINT_DATA" << " " << nx * ny * 1 << endl;
	outfile << "SCALARS" << " " << scalar_attribute_1 << " " << "double 1" << endl;
	outfile << "LOOKUP_TABLE default" << endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			outfile << Phi[i][j] << endl;
		}
	}
	outfile << "SCALARS" << " " << scalar_attribute_2 << " " << "double 1" << endl;
	outfile << "LOOKUP_TABLE default" << endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			outfile << (area[i][j] == interface ? 1.0 : 0.0) << endl;
		}
	}
	outfile << "VECTORS" << " " << vector_attribute_1 << " " << "float" << endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			outfile << u_vector_x[i][j] << " " << v_vector_y[i][j] << " " << 0 << endl;
		}
	}
	outfile.close();
	return;
}

void bounceBack(
	const Region::regiontype(&region)[nx][ny],
	double(&f)[nx][ny][9],
	double(&u)[nx][ny],
	double(&v)[nx][ny]) {
	//����Ϊ bounceBack(Region::region,f,u,v);
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			if (region[i][j] == Region::FB) {

				if (region[i - 1][j] == Region::SB) {
					f[i][j][1] = f[i - 1][j][3];
				}
				if (region[i][j - 1] == Region::SB) {
					f[i][j][2] = f[i][j - 1][4];
				}
				if (region[i + 1][j] == Region::SB) {
					f[i][j][3] = f[i + 1][j][1];
				}
				if (region[i][j + 1] == Region::SB) {
					f[i][j][4] = f[i][j + 1][2];
				}
				if (region[i - 1][j - 1] == Region::SB) {
					f[i][j][5] = f[i - 1][j - 1][7];
				}
				if (region[i + 1][j - 1] == Region::SB) {
					f[i][j][6] = f[i + 1][j - 1][8];
				}
				if (region[i + 1][j + 1] == Region::SB) {
					f[i][j][7] = f[i + 1][j + 1][5];
				}
				if (region[i - 1][j + 1] == Region::SB) {
					f[i][j][8] = f[i - 1][j + 1][6];
				}

			}

		}
	}
	return;
}

void diffuse(double(&Phi)[nx][ny]) {
	//�ೡ����ķ�ɢ��
	for (int iter = 1; iter <= 3; iter++) {
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
					Phi[i][j] = 4.0 / 9.0 * Phi[i][j] + 1.0 / 9.0 * (Phi[i + 1][j] + Phi[i][j + 1] + Phi[i - 1][j] + Phi[i][j - 1])
						+ 1.0 / 36.0 * (Phi[i + 1][j + 1] + Phi[i + 1][j - 1] + Phi[i - 1][j + 1] + Phi[i - 1][j - 1]);
				}
			}
		}
	}
	return;
}

void cal_Phi(double(&rho_r)[nx][ny], double(&rho_b)[nx][ny], double(&rho)[nx][ny]) {
	//������ɫ��������������
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {

				Phi[i][j] = (rho_r[i][j] - rho_b[i][j]) / (rho_r[i][j] + rho_b[i][j]);

			}
		}
	}
	return;
}

void calculate_temp_moment_eq(double& rho, double& u, double& v) {
	//������ʱƽ���
	//�β�rho��Ӧʵ��rho[i][j]
	//�β�u��Ӧʵ��u[i][j]
	//�β�v��Ӧʵ��v[i][j]
	temp_moment_eq[0] = rho;
	temp_moment_eq[1] = -2.0 * rho + 3.0 * rho * (u * u + v * v);
	temp_moment_eq[2] = -3.0 * rho * (u * u + v * v) + rho;
	//temp_moment_eq[2] = 9.0 * rho * u * u * v * v - 3.0 * rho * (u * u + v * v) + rho;
	temp_moment_eq[3] = rho * u;
	temp_moment_eq[4] = -rho * u;
	//temp_moment_eq[4] = rho * u * (3.0 * v * v - 1);
	temp_moment_eq[5] = rho * v;
	temp_moment_eq[6] = -rho * v;
	//temp_moment_eq[6] = rho * v * (3.0 * u * u - 1);
	temp_moment_eq[7] = rho * (u * u - v * v);
	temp_moment_eq[8] = rho * u * v;
	return;
}

void collision_single_phase(
	double(&f)[nx][ny][Q],
	const double(&M)[Q][Q],
	const double(&InvM)[Q][Q],
	double(&S)[Q],
	double(&temp_moment)[Q],
	double(&temp_moment_post_collison)[Q],
	double(&rho)[nx][ny],
	double(&u)[nx][ny],
	double(&v)[nx][ny]
) {
	// �������ӿռ�ת�����ؿռ䡢�ؿռ����ײ���ؿռ�ת�������ӿռ�
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {

				// ת�����ؿռ�
				for (int k = 0; k < 9; k++) {
					temp_moment[k] = M[k][0] * f[i][j][0] + M[k][1] * f[i][j][1] + M[k][2] * f[i][j][2] +
						M[k][3] * f[i][j][3] + M[k][4] * f[i][j][4] + M[k][5] * f[i][j][5] +
						M[k][6] * f[i][j][6] + M[k][7] * f[i][j][7] + M[k][8] * f[i][j][8];
				}

				//������ʱƽ���
				calculate_temp_moment_eq(rho[i][j], u[i][j], v[i][j]);//�����rhoΪ�����Ĳ���������ȫ�ֱ���

				// �ؿռ����ײ
				//S[7]��S[8]�ڽ��洦�Ǳ仯��
				for (int k = 0; k < 7; k++) {
					temp_moment_post_collison[k] = (1 - S[k]) * temp_moment[k] + S[k] * temp_moment_eq[k];
					//temp_moment_post_collison[k] = (1 - 0.8 * omega(i, j)) * temp_moment[k] + 0.8 * omega(i, j) * temp_moment_eq[k];
				}
				temp_moment_post_collison[7] = (1 - omega(i, j)) * temp_moment[7] + omega(i, j) * temp_moment_eq[7];
				temp_moment_post_collison[8] = (1 - omega(i, j)) * temp_moment[8] + omega(i, j) * temp_moment_eq[8];

				// ת�������ӿռ�,�ֲ����� f �����˾ؿռ��ڵ���ײ
				for (int k = 0; k < 9; k++) {
					f[i][j][k] = InvM[k][0] * temp_moment_post_collison[0] + InvM[k][1] * temp_moment_post_collison[1] +
						InvM[k][2] * temp_moment_post_collison[2] + InvM[k][3] * temp_moment_post_collison[3] +
						InvM[k][4] * temp_moment_post_collison[4] + InvM[k][5] * temp_moment_post_collison[5] +
						InvM[k][6] * temp_moment_post_collison[6] + InvM[k][7] * temp_moment_post_collison[7] +
						InvM[k][8] * temp_moment_post_collison[8];
				}

			}
		}
	}
	return;
}

void calculate_temp_moment_F(double(&F)[2], double& u, double& v) {
	//�β�F[2]��Ӧʵ��F[i][j]
	//�β�u��Ӧʵ��u[i][j]
	//�β�v��Ӧʵ��v[i][j]
	//calculate_temp_moment_F(F[i][j],u[i][j],v[i][j]);
	temp_moment_F[0] = 0;
	temp_moment_F[1] = 6.0 * (F[0] * u + F[1] * v);
	temp_moment_F[2] = -6.0 * (F[0] * u + F[1] * v);
	temp_moment_F[3] = F[0];
	temp_moment_F[4] = -F[0];
	temp_moment_F[5] = F[1];
	temp_moment_F[6] = -F[1];
	temp_moment_F[7] = 2.0 * (F[0] * u - F[1] * v);
	temp_moment_F[8] = F[1] * u + F[0] * v;
	return;
}

void sourceTermPrecondition() {
	//ʹ�ø�Դ�������ʩ�ӱ������� (Ԥ����)

	//��һ�������ೡ�ݶȵ�λ����
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
				if (area[i][j] == interface) {
					n[i][j][0] = F_x[i][j] / (F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j]);
					n[i][j][1] = F_y[i][j] / (F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j]);
				}
			}
		}
	}

	//�ڶ������� ���� ������ ���������� ���Ĳ� Դ�� MF
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
				if (area[i][j] == interface) {
					//����
					kappa[i][j] = -(1.0 / 12.0 * (n[i + 1][j + 1][0] - n[i - 1][j + 1][0] + n[i + 1][j - 1][0] - n[i - 1][j - 1][0] +
						n[i + 1][j + 1][1] - n[i + 1][j - 1][1] + n[i - 1][j + 1][1] - n[i - 1][j - 1][1]) +
						1.0 / 3.0 * (n[i + 1][j][0] - n[i - 1][j][0] + n[i][j + 1][1] - n[i][j - 1][1]));

					//����������
					F[i][j][0] = 0.5 * theta * kappa[i][j] * F_x[i][j];
					F[i][j][1] = 0.5 * theta * kappa[i][j] * F_y[i][j];
					//Դ�� M^(-1) * (I  0.5 * S) * MF 
					//����main�м���
				}
				else
				{
					//ע�⽫���ڽ���ı�����������Ϊ0
					F[i][j][0] = 0.0;
					F[i][j][1] = 0.0;
				}
			}
		}
	}

	return;
}

void debug(int kk, double(&N_r)[nx][ny][Q], double(&N_b)[nx][ny][Q], int step) {
	//debug(1,N_r,N_b,step);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < Q; k++) {
				if (isnan(N_r[i][j][k]) != 0 || isnan(N_b[i][j][k]) != 0 || isnan(Phi[i][j]) != 0) {
					cout << step << "	" << kk << endl;
					return;
				}
			}
		}
	}
}

int main() {
	Region::setRegion::setProcess();//����region
	Region::cal_ns();//�������߽�ĵ�λ��������ָ�����
	inoutBoundary::caculate_uInletDistribution();//������ڵ��ٶȷֲ�
	//��ʼ��
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] == Region::SL || Region::region[i][j] == Region::SB) {
				area[i][j] = other;
			}
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
				area[i][j] = blue;
			}
		}
	}

	//���Һ�ε�λ�ò���
	//double length;
	//double obstacle_cylinder_centerx[2] = { 3 * nx / 7.0,4 * nx / 7.0 };
	//double obstacle_cylinder_certery[2] = { 0,0 };
	//double obstacle_cylinder_r = (nx < ny ? nx : ny) / 5.0;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
				//for (int k = 0; k < 2; k++) {
				//	length = sqrt(double(pow(i - obstacle_cylinder_centerx[k], 2) + pow(j - obstacle_cylinder_certery[k], 2)));
				//	if (length < obstacle_cylinder_r)
				//	{
				//		area[i][j] = red;
				//	} 
				//}
				/*ppt��ʹ�õ��������ʼ�ֲ�!!!!!!!!!!!!!!!!!
				if (i > 2.7 * nx / 7.0 && i < 4.3 * nx / 7.0 && j < ny / 3.0) {
					area[i][j] = red;
				}
				*/
				/*if (i > 3 * nx / 7.0 && i < 4 * nx / 7.0 && j < ny / 4.0) {
					area[i][j] = red;
				}*/
				if (i < 10) {
					area[i][j] = red;
				}
				//if (i < 1.0 * nx / 7.0 && j < ny / 4.0) {
				//	area[i][j] = red;
				//}


			}
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (area[i][j] == red)
			{
				rho_r[i][j] = rho_r0;
				rho_b[i][j] = 0;
			}
			if (area[i][j] == blue)
			{
				rho_r[i][j] = 0;
				rho_b[i][j] = rho_b0;
			}
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
				rho[i][j] = rho_r[i][j] + rho_b[i][j];
			}
		}
	}

	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			for (int k = 0; k < Q; k++) {
				N_r[i][j][k] = w[k] * rho_r[i][j];//��ʼ��
				N_b[i][j][k] = w[k] * rho_b[i][j];
			}
		}
	}
	double temp = 0;
	//��ʼ
	for (int step = 0; step < 8001; step++) {
		if (step % 10 == 0) {
			temp = 0;
			//cout << "��" << step << "��" << endl;
			cout << "��" << step << "��";
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					uu = u[i][j] * u[i][j] + v[i][j] * v[i][j];
					if (temp < sqrt(uu)) {
						temp = sqrt(uu);
					}
				}
			}
			cout << "����ٶ�Ϊ:" << temp << endl;
		}
		//������ײ����
		cal_Phi(rho_r, rho_b, rho);//������ɫ��,��׽����
		// �������ӿռ�ת�����ؿռ䡢�ؿռ����ײ���ؿռ�ת�������ӿռ�

		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				for (int k = 0; k < Q; k++) {
					N[i][j][k] = N_r[i][j][k] + N_b[i][j][k];
				}
			}
		}

		collision_single_phase(N, M, InvM, S, temp_moment, temp_moment_post_collison, rho, u, v);

		//����1��ʼ==================================================
		debug(1, N_r, N_b, step);
		//����1����==================================================

		//�Ŷ�����
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				if (Region::region[i][j] == Region::FL) {

					//������ɫ�ݶ�
					F_x[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i - 1][j + 1])
						+ 2.0 / 3.0 * (Phi[i + 1][j] - Phi[i - 1][j]) + 1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i - 1][j - 1]));
					F_y[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i + 1][j - 1])
						+ 2.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i - 1][j + 1] - Phi[i - 1][j - 1]));

					//��������1 ���ೡ�ݶȵĴ�С����aqrt(delta2)ʱ��Ϊ����
					if (F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j] > delta2) {
						area[i][j] = interface;
					}
					else {
						area[i][j] = other;
					}
				}
			}
		}

		//ʪ��߽�����

		wettingBoundary::wetting_boundary(F_x, F_y);

		//����2��ʼ==================================================
		debug(2, N_r, N_b, step);
		//����2����==================================================

		//���������������Դ���Ԥ����
		sourceTermPrecondition();

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {

				if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
					if (area[i][j] == interface) {

						calculate_temp_moment_F(F[i][j], u[i][j], v[i][j]);

						for (int k = 0; k < 7; k++) {
							temp_moment_post_collison_F[k] = (1 - S[k]) * temp_moment_F[k];
							//temp_moment_post_collison_F[k] = (1 - 0.8 * omega(i, j)) * temp_moment_F[k];
						}
						temp_moment_post_collison_F[7] = (1 - omega(i, j)) * temp_moment_F[7];
						temp_moment_post_collison_F[8] = (1 - omega(i, j)) * temp_moment_F[8];

						// ת�������ӿռ�,�ֲ����� f �����˾ؿռ��ڵ���ײ
						for (int k = 0; k < 9; k++) {

							N[i][j][k] += InvM[k][0] * temp_moment_post_collison_F[0] + InvM[k][1] * temp_moment_post_collison_F[1] +
								InvM[k][2] * temp_moment_post_collison_F[2] + InvM[k][3] * temp_moment_post_collison_F[3] +
								InvM[k][4] * temp_moment_post_collison_F[4] + InvM[k][5] * temp_moment_post_collison_F[5] +
								InvM[k][6] * temp_moment_post_collison_F[6] + InvM[k][7] * temp_moment_post_collison_F[7] +
								InvM[k][8] * temp_moment_post_collison_F[8];

						}

					}
				}

			}
		}

		//����3��ʼ==================================================
		debug(3, N_r, N_b, step);
		//����3����==================================================

		//����ɫ����

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
					if (area[i][j] == interface) {
						//һ����˵������Ӧ�����һ���жϣ�����interface��ʹ��
						FF = F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j];
						for (int k = 0; k < Q; k++) {
							Fc = F_x[i][j] * cx[k] + F_y[i][j] * cy[k];

							N_r[i][j][k] = rho_r[i][j] / rho[i][j] * N[i][j][k] +
								beta * rho_r[i][j] * rho_b[i][j] / rho[i][j] * w[k] * cosphiImprove(k, Fc, FF);
							N_b[i][j][k] = rho_b[i][j] / rho[i][j] * N[i][j][k] -
								beta * rho_r[i][j] * rho_b[i][j] / rho[i][j] * w[k] * cosphiImprove(k, Fc, FF);
						}
					}
				}
			}
		}

		//����4��ʼ==================================================
		debug(4, N_r, N_b, step);
		//����4����==================================================

		//��
		streaming(N_r);
		streaming(N_b);

		//����5��ʼ==================================================
		debug(5, N_r, N_b, step);
		//����5����==================================================

		//�߽�
		bounceBack(Region::region, N_r, u, v);
		bounceBack(Region::region, N_b, u, v);

		//ʹ���ֲܷ�������������ڳ��ڱ߽�
		inoutBoundary::inletVelocity(N, N_r, N_b, rho, rho_r, rho_b, u, v, uInletDistribution);
		inoutBoundary::outletPressure(N, N_r, N_b, rho, rho_r, rho_b, u, v, rhoOut);



		//����6��ʼ==================================================
		debug(6, N_r, N_b, step);
		//����6����==================================================

		//�ܶ�
		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
					rho_r[i][j] = N_r[i][j][0] + N_r[i][j][1] + N_r[i][j][2] + N_r[i][j][3] +
						N_r[i][j][4] + N_r[i][j][5] + N_r[i][j][6] + N_r[i][j][7] + N_r[i][j][8];
					rho_b[i][j] = N_b[i][j][0] + N_b[i][j][1] + N_b[i][j][2] + N_b[i][j][3] +
						N_b[i][j][4] + N_b[i][j][5] + N_b[i][j][6] + N_b[i][j][7] + N_b[i][j][8];
					rho[i][j] = rho_r[i][j] + rho_b[i][j];
				}
			}
		}

		//�ٶ�
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
					u[i][j] = ((N_r[i][j][1] + N_r[i][j][5] + N_r[i][j][8]) - (N_r[i][j][3] + N_r[i][j][6] + N_r[i][j][7]) +
						(N_b[i][j][1] + N_b[i][j][5] + N_b[i][j][8]) - (N_b[i][j][3] + N_b[i][j][6] + N_b[i][j][7]) + 0.5 * F[i][j][0]) / rho[i][j];
					v[i][j] = ((N_r[i][j][2] + N_r[i][j][5] + N_r[i][j][6]) - (N_r[i][j][4] + N_r[i][j][7] + N_r[i][j][8]) +
						(N_b[i][j][2] + N_b[i][j][5] + N_b[i][j][6]) - (N_b[i][j][4] + N_b[i][j][7] + N_b[i][j][8]) + 0.5 * F[i][j][1]) / rho[i][j];
				}
			}
		}


		//if (step % interval == 0) {
		//	out_file("Phi" + to_string(step / interval) + ".vtk", "Phi", "velocity", u, v, Phi);
		//}
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (Region::region[i][j] == Region::SB || Region::region[i][j] == Region::SL) {
					Phi[i][j] = 0;
				}
			}
		}
		if (step % interval == 0) {
			out_file("new" + to_string(step / interval) + ".vtk", "Phi", "area", "velocity", u, v, Phi, area);
		}
		//����7��ʼ==================================================
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int k = 0; k < Q; k++) {
					if (isnan(N_r[i][j][k]) != 0 || isnan(N_b[i][j][k]) != 0) {
						cout << step << "	7" << endl;
						return -7;
					}
				}
			}
		}
		//����7����==================================================
	}
	return 0;
}