#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <fstream>

//����LBMCGM-2020-11-2�Ľ���ϣ��ʹ��ѹ���߽���������֤���ȶ�ǰ�ص��γ�
//������ļ�������ѹ��,��rho
//��ʹ�õı߽���������"Lattice Boltzmann simulation of pressure-driven two-phase flows in capillary tube and porous medium"

using namespace std;
const double M[9][9] = {
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1},
	{-4,-1,-1,-1,-1, 2, 2, 2, 2},
	{ 4,-2,-2,-2,-2, 1, 1, 1, 1},
	{ 0, 1, 0,-1, 0, 1,-1,-1, 1},
	{ 0,-2, 0, 2, 0, 1,-1,-1, 1},
	{ 0, 0, 1, 0,-1, 1, 1,-1,-1},
	{ 0, 0,-2, 0, 2, 1, 1,-1,-1},
	{ 0, 1,-1, 1,-1, 0, 0, 0, 0},
	{ 0, 0, 0, 0, 0, 1,-1, 1,-1} };

const double InvM[9][9] = {
	{1.0 / 9.0, -1.0 / 9.0,  1.0 / 9.0,    0.0,    0.0,    0.0,   0.0,   0.0,  0.0},
	{1.0 / 9.0, -1.0 / 36.0,-1.0 / 18.0,1.0 / 6.0,-1.0 / 6.0,0.0,0.0,1.0 / 4.0,0.0},
	{1.0 / 9.0, -1.0 / 36.0,-1.0 / 18.0,0.0,0.0,1.0 / 6.0,-1.0 / 6.0,-1.0 / 4.0,0.0},
	{1.0 / 9.0, -1.0 / 36.0,-1.0 / 18.0,-1.0 / 6.0,1.0 / 6.0,0.0,0.0,1.0 / 4.0,0.0},
	{1.0 / 9.0, -1.0 / 36.0,-1.0 / 18.0,0.0,0.0,-1.0 / 6.0,1.0 / 6.0,-1.0 / 4.0,0.0},
	{1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0,1.0 / 6.0,1.0 / 12.0,1.0 / 6.0,1.0 / 12.0,0.0,1.0 / 4.0},
	{1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0,-1.0 / 6.0,-1.0 / 12.0,1.0 / 6.0,1.0 / 12.0,0.0,-1.0 / 4.0},
	{1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0,-1.0 / 6.0,-1.0 / 12.0,-1.0 / 6.0,-1.0 / 12.0,0.0,1.0 / 4.0},
	{1.0 / 9.0,  1.0 / 18.0, 1.0 / 36.0,1.0 / 6.0,1.0 / 12.0,-1.0 / 6.0,-1.0 / 12.0,0.0,-1.0 / 4.0} };
const int nx = 100;
const int ny = 50;
const int Q = 9;
double N_r[nx][ny][Q]{ 0 };
double N_b[nx][ny][Q]{ 0 };
double N[nx][ny][Q]{ 0 };
//double Neq_r[nx][ny][Q]{ 0 };
//double Neq_b[nx][ny][Q]{ 0 };
//double Neq[nx][ny][Q]{ 0 };
double rho_r[nx][ny]{ 0 };
double rho_b[nx][ny]{ 0 };
double rho[nx][ny]{ 0 };
double u[nx][ny]{ 0 };
double v[nx][ny]{ 0 };
double rhoInlet = 1.01;//����ܶȣ�����ѹ����ڱ߽�
double rhoOutlet = 0.99;//�����ܶȣ�����ѹ�����ڱ߽�
double uInlet = 0.05;//�������
double uInletDistribution[ny];
double vInlet = 0.0;
double Phi[nx][ny]{ 0 };//�ೡ
double F_x[nx][ny]{ 0 };//�ೡ�ݶ�
double F_y[nx][ny]{ 0 };
double cx[9]{ 0,1,0,-1,0,1,-1,-1,1 };
double cy[9]{ 0,0,1,0,-1,1,1,-1,-1 };
double w[9] = { 4.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0 };
double rho_r0 = 1;//��ʼ�ܶ�
double rho_b0 = 1;
double nu_r = 0.4;//ճ��ϵ����ճ��Խ�󣬳�ԥ�ʾ�Խ��,����ƽ����ٶ�Ҳ��Խ��,���Ҳ��Ӱ��α���������������ճ�Զ����ʱ��α�������С
//nu_r�̶����䣬����nu_b�����ճ�ͱ�,�������nu_b[1.0/6000.0,1000.0/6.0]֮�䣬������ճ��ϵ����Ϊ 0.001~1000
double nu_b = 0.4;
double contactAngle = 135.0 / 180.0 * 3.141592653;//�Ӵ���!!!!!!!!!!!!!!!!!
double A = 0.01;
double beta = 1.0;
double delta2 = 1e-4;//�����ೡ�ݶȿ��ƽ�����,Ŀǰ���ֵ��ȡ���Ǹ���ͼ���Գ�����1e-4
double uu, cu, FF, Fc;
// S�ԽǾ��� �ɳ�ϵ������ s7��s8��BGK�е�omega��ͬ //���Ե���һ�µ�ϵ��������ȶ�
double s0 = 1.0, s1 = 1.64, s2 = 1.54, s3 = 1.0, s4 = 1.9, s5 = 1.0, s6 = 1.9, s7, s8;//s7 s8��仯//s4��s6Ҳ����Ҫ�Ĳ��������ڿ����ȶ� huang���õ���1.2 liu���õ���1.9
double S[9] = { s0,s1,s2,s3,s4,s5,s6,s7,s8 };
// ƽ��طֲ�����
double temp_moment[9]{ 0 };
double temp_moment_eq[9]{ 0 };
double temp_moment_post_collison[9]{ 0 };

enum areatype { interface, red, blue, contactline, solid };
areatype area[nx][ny];

int interval = 20;

namespace Region {
	enum regiontype { F, FB, FL, S, SB, SL, inlet, outlet, FB_inlet, FB_outlet };//һ����FB,FL,SB,SL, inlet, outlet,FB_inlet,FB_outlet 8����ɲ���
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
		else if (region[i][j] == FB || region[i][j] == FL || region[i][j] == FB_inlet || region[i][j] == FB_outlet) {
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
				if (region[i][j] == FB || region[i][j] == FB_inlet || region[i][j] == FB_outlet) {
					ns[i][j][0] = tempx(i, j) / sqrt(tempx(i, j) * tempx(i, j) + tempy(i, j) * tempy(i, j));
					ns[i][j][1] = tempy(i, j) / sqrt(tempx(i, j) * tempx(i, j) + tempy(i, j) * tempy(i, j));
				}
			}
		}
		//����ns
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (region[i][j] != FB && region[i][j] != FB_inlet && region[i][j] != FB_outlet) {
					ns[i][j][0] = 0;
					ns[i][j][1] = 0;
				}
			}
		}
		//��ư����±߽��ns
		for (int i = 0; i < nx; i++) {
			ns[i][2][0] = 0;
			ns[i][2][1] = -1;
			ns[i][ny - 3][0] = 0;
			ns[i][ny - 3][1] = 1;
		}
		return;
	}

}

namespace Region {

	//����Region::setProcess();

	//���2*4��Բ���뾶Ϊ6
	const int x = 4;
	const int y = 2;
	double r = 6;
	double length;
	long long site_x[x][y];
	long long site_y[x][y];
	void setS() {

		for (int i = 0; i < x; i++) {
			for (int j = 0; j < y; j++) {
				site_x[i][j] = nx / double(x + 1) * (double(i) + 1.0);
				site_y[i][j] = ny / double(y + 1) * (double(j) + 1.0);
			}
		}		

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (j == 0 || j == 1 || j == ny - 1 || j == ny - 2) {
					region[i][j] = S;
				}

				for (int ii = 0; ii < x; ii++) {
					for (int jj = 0; jj < y; jj++) {
						length = sqrt(pow(i - site_x[ii][jj], 2) + pow(j - site_y[ii][jj], 2));
						if (length < r) {
							region[i][j] = S;
						}
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
				if (region[i][j] == S) {
					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (region[i + x][j + y] == F) {
								region[i + x][j + y] = FB;
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
				if (region[i][j] == FB) {
					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (region[i + x][j + y] == S) {
								region[i + x][j + y] = SB;
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

	void setinlet_outlet() {
		for (int j = 2; j <= ny - 3; j++) {
			region[0][j] = inlet;
			region[nx - 1][j] = outlet;
		}
		return;
	}

	void setFBinlet_FBoutlet() {
		region[0][2] = FB_inlet;
		region[0][ny - 3] = FB_inlet;
		region[nx - 1][2] = outlet;
		region[nx - 1][ny - 3] = FB_outlet;
		return;
	}

	void setProcess() {
		setS();
		setF();
		setFBFL();
		setSBSL();
		setinlet_outlet();
		setFBinlet_FBoutlet();
		return;
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
					Phi[i][j] = 0;
					temp1 = 0;
					temp2 = 0;
					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (Region::region[i + x][j + y] == Region::FB || Region::region[i + x][j + y] == Region::FB_inlet || Region::region[i + x][j + y] == Region::FB_outlet) {
								temp1 += w(x * x + y * y) * Phi[i + x][j + y];
								temp2 += w(x * x + y * y);
							}
						}
					}
					for (int x = -1; x <= 1; x++) {
						for (int y = -1; y <= 1; y++) {
							if (Region::region[i + x][j + y] == Region::FB || Region::region[i + x][j + y] == Region::FB_inlet || Region::region[i + x][j + y] == Region::FB_outlet) {
								Phi[i][j] = temp1 / temp2;
							}
						}
					}
				}
			}
		}

		//���ұ߽�ǵ㴦S���ೡ
		Phi[0][1] = (w(1) * Phi[0][2] + w(2) * Phi[1][2]) / (w(1) + w(2));
		Phi[0][ny - 2] = (w(1) * Phi[0][ny - 3] + w(2) * Phi[1][ny - 3]) / (w(1) + w(2));
		Phi[nx - 1][1] = (w(1) * Phi[nx - 1][2] + w(2) * Phi[nx - 2][2]) / (w(1) + w(2));
		Phi[nx - 1][ny - 2] = (w(1) * Phi[nx - 1][ny - 3] + w(2) * Phi[nx - 2][ny - 3]) / (w(1) + w(2));

		double n1x, n1y, n2x, n2y, Fxn, Fyn, D1, D2;
		//�ڶ������õ�һ�������SB�����ೡPhi����FB�����ೡ�ݶ�F_x,F_y
		for (int i = 0; i < nx; i++) {
			for (int j = 1; j < ny - 1; j++) {

				if (Region::region[i][j] == Region::FB || Region::region[i][j] == Region::FB_inlet || Region::region[i][j] == Region::FB_outlet) {

					if (i != 0 && i != nx - 1) {
						F_x[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i - 1][j + 1])
							+ 2.0 / 3.0 * (Phi[i + 1][j] - Phi[i - 1][j]) + 1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i - 1][j - 1]));

						F_y[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i + 1][j - 1])
							+ 2.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i - 1][j + 1] - Phi[i - 1][j - 1]));
					}
					else if (i == 0 && j == 2) {
						F_x[i][j] = 2.0 / 3.0 * (Phi[i + 1][j] - Phi[i][j]) + 1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i][j + 1]) +
							1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i][j - 1]);

						F_y[i][j] = 1.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i + 1][j - 1]);
					}
					else if (i == 0 && j == ny - 3) {
						F_x[i][j] = 2.0 / 3.0 * (Phi[i + 1][j] - Phi[i][j]) + 1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i][j + 1]) +
							1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i][j - 1]);

						F_y[i][j] = 1.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i + 1][j - 1]);
					}
					else if (i == nx - 1 && j == 2) {
						F_x[i][j] = 2.0 / 3.0 * (Phi[i][j] - Phi[i - 1][j]) + 1.0 / 6.0 * (Phi[i][j + 1] - Phi[i - 1][j + 1]) +
							1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i][j - 1]);

						F_y[i][j] = 1.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i - 1][j + 1] - Phi[i - 1][j - 1]);
					}
					else if (i == nx - 1 && j == ny - 3) {
						F_x[i][j] = 2.0 / 3.0 * (Phi[i][j] - Phi[i - 1][j]) + 1.0 / 6.0 * (Phi[i][j + 1] - Phi[i - 1][j + 1]) +
							1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i][j - 1]);

						F_y[i][j] = 1.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i - 1][j + 1] - Phi[i - 1][j - 1]);
					}

					//��������2 ���ೡ�ݶȵĴ�С����aqrt(delta2)ʱ��Ϊ����
					if (F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j] > delta2) {
						area[i][j] = contactline;

						//������������FB������ೡ�ݶ�F_x,F_y��С���䣬���ݽӴ�����FB���Ḻ́ڷ���ns��������

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
					else if (Phi[i][j] > 0) {
						area[i][j] = red;
					}
					else
						area[i][j] = blue;
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

	void inletVelocity() {

		double temp_rho = 0;
		int count = 0;

		//�ϲ��ֲ�����
		for (int j = 0; j < ny; j++) {
			if (Region::region[0][j] == Region::inlet || Region::region[0][j] == Region::FB_inlet)
				for (int k = 0; k < Q; k++) {
					N[0][j][k] = N_r[0][j][k] + N_b[0][j][k];
				}
		}

		//��ڵĵ�λ����Ϊ(1,0)ʱ���ٶȱ߽磬(1,0)ָ������
		for (int j = 0; j < ny; j++) {
			if (Region::region[0][j] == Region::inlet) {
				rho[0][j] = (N[0][j][0] + N[0][j][2] + N[0][j][4] + 2.0 * (N[0][j][3] + N[0][j][6] + N[0][j][7])) / (1 - uInletDistribution[j]);
				N[0][j][1] = N[0][j][3] + 2.0 * rho[0][j] * uInletDistribution[j] / 3.0;
				N[0][j][5] = N[0][j][7] - 0.5 * (N[0][j][2] - N[0][j][4]) + rho[0][j] * uInletDistribution[j] / 6.0;
				N[0][j][8] = N[0][j][6] + 0.5 * (N[0][j][2] - N[0][j][4]) + rho[0][j] * uInletDistribution[j] / 6.0;
				temp_rho += rho[0][j];
				count++;
			}
		}

		temp_rho = temp_rho / count;

		//FB_inlet��Ҫ�ر���

		N[0][2][1] = N[0][2][3]; N[0][2][2] = N[0][2][4]; N[0][2][5] = N[0][2][7];
		N[0][2][6] = N[0][2][8] = 0.5 * (temp_rho - (N[0][2][0] + N[0][2][1] + N[0][2][2] + N[0][2][3] + N[0][2][4] + N[0][2][5] + N[0][2][7]));

		N[0][ny - 3][1] = N[0][ny - 3][3]; N[0][ny - 3][4] = N[0][ny - 3][2]; N[0][ny - 3][8] = N[0][ny - 3][6];
		N[0][ny - 3][5] = N[0][ny - 3][7] = 0.5 * (temp_rho - (N[0][ny - 3][0] + N[0][ny - 3][1] + N[0][ny - 3][2] + N[0][ny - 3][3] + N[0][ny - 3][4] + N[0][ny - 3][6] + N[0][ny - 3][8]));

		//�����ܶȷֲ��ֽ�ֲ�����
		for (int j = 0; j < ny; j++) {
			if (Region::region[0][j] == Region::inlet || Region::region[0][j] == Region::FB_inlet)
				for (int k = 0; k < Q; k++) {
					N_r[0][j][k] = N[0][j][k] * rho_r[0][j] / (rho_r[0][j] + rho_b[0][j]);/////=================-------------===============
					N_b[0][j][k] = N[0][j][k] - N_r[0][j][k];
				}
		}

		return;
	}

	void ouletNeumann(double(&f)[nx][ny][Q]) {//------------=======================------------------==================-

		for (int j = 0; j < ny; j++) {
			if (Region::region[nx - 1][j] == Region::outlet) {
				for (int q = 0; q < Q; q++) {
					f[nx - 1][j][q] = f[nx - 2][j][q];
				}
			}
			else if (j == 2) {
				for (int q = 0; q < Q; q++) {
					f[nx - 1][j][q] = f[nx - 2][j][q];
				}
			}
			else if (j == ny - 3) {
				for (int q = 0; q < Q; q++) {
					f[nx - 1][j][q] = f[nx - 2][j][q];
				}
			}
		}

		return;
	}

	//Ϊ�˳��������Ե�Ĳ��ȶ��ԣ�����ѹ���߽�����
	void inletPressure() {

		double rhoinlet_ux;

		//�ϲ��ֲ�����
		for (int j = 0; j < ny; j++) {
			if (Region::region[0][j] == Region::inlet || Region::region[0][j] == Region::FB_inlet)
				for (int k = 0; k < Q; k++) {
					N[0][j][k] = N_r[0][j][k] + N_b[0][j][k];
				}
		}

		//��ڵĵ�λ����Ϊ(1,0)ʱ���ٶȱ߽磬(1,0)ָ������
		for (int j = 0; j < ny; j++) {
			if (Region::region[0][j] == Region::inlet) {
				rhoinlet_ux = rhoInlet - (N[0][j][0] + N[0][j][2] + N[0][j][4]) - 2 * (N[0][j][3] + N[0][j][6] + N[0][j][7]);
				N[0][j][1] = N[0][j][3] + 2.0 / 3.0 * rhoinlet_ux;
				N[0][j][5] = N[0][j][7] - 0.5 * (N[0][j][2] - N[0][j][4]) + rhoinlet_ux / 6.0;
				N[0][j][8] = N[0][j][6] + 0.5 * (N[0][j][2] - N[0][j][4]) + rhoinlet_ux / 6.0;
			}
		}

		//FB_inlet��Ҫ�ر���

		N[0][2][1] = N[0][2][3]; N[0][2][2] = N[0][2][4]; N[0][2][5] = N[0][2][7];
		N[0][2][6] = N[0][2][8] = 0.5 * (rhoInlet - (N[0][2][0] + N[0][2][1] + N[0][2][2] + N[0][2][3] + N[0][2][4] + N[0][2][5] + N[0][2][7]));

		N[0][ny - 3][1] = N[0][ny - 3][3]; N[0][ny - 3][4] = N[0][ny - 3][2]; N[0][ny - 3][8] = N[0][ny - 3][6];
		N[0][ny - 3][5] = N[0][ny - 3][7] = 0.5 * (rhoInlet - (N[0][ny - 3][0] + N[0][ny - 3][1] + N[0][ny - 3][2] + N[0][ny - 3][3] + N[0][ny - 3][4] + N[0][ny - 3][6] + N[0][ny - 3][8]));

		//�����ܶȷֲ��ֽ�ֲ�����
		for (int j = 0; j < ny; j++) {
			if (Region::region[0][j] == Region::inlet || Region::region[0][j] == Region::FB_inlet)
				for (int k = 0; k < Q; k++) {
					N_r[0][j][k] = N[0][j][k] * rho_r[0][j] / (rho_r[0][j] + rho_b[0][j]);/////=================-------------===============
					N_b[0][j][k] = N[0][j][k] - N_r[0][j][k];
				}
		}


		return;
	}

	void outletPressure() {

		double rhooutlet_ux;

		//�ϲ��ֲ�����
		for (int j = 0; j < ny; j++) {
			if (Region::region[nx - 1][j] == Region::outlet || Region::region[nx - 1][j] == Region::FB_outlet)
				for (int k = 0; k < Q; k++) {
					N[nx - 1][j][k] = N_r[nx - 1][j][k] + N_b[nx - 1][j][k];
				}
		}

		//���ڵĵ�λ����Ϊ(1,0)ʱ���ٶȱ߽磬(1,0)ָ������
		for (int j = 0; j < ny; j++) {
			if (Region::region[nx - 1][j] == Region::outlet) {
				rhooutlet_ux = (N[nx - 1][j][0] + N[nx - 1][j][2] + N[nx - 1][j][4]) + 2 * (N[nx - 1][j][1] + N[nx - 1][j][5] + N[nx - 1][j][8]) - rhoOutlet;
				N[nx - 1][j][3] = N[nx - 1][j][1] - 2.0 / 3.0 * rhooutlet_ux;
				N[nx - 1][j][6] = N[nx - 1][j][8] - rhooutlet_ux / 6.0 - 0.5 * (N[nx - 1][j][2] - N[nx - 1][j][4]);
				N[nx - 1][j][7] = N[nx - 1][j][5] - rhooutlet_ux / 6.0 + 0.5 * (N[nx - 1][j][2] - N[nx - 1][j][4]);
			}
		}

		//FB_outlet��Ҫ�ر���

		N[nx - 1][2][2] = N[nx - 1][2][4]; N[nx - 1][2][3] = N[nx - 1][2][1]; N[nx - 1][2][6] = N[nx - 1][2][8];
		N[nx - 1][2][5] = N[nx - 1][2][7] = 0.5 * (rhoOutlet - (N[nx - 1][2][0] + N[nx - 1][2][1] + N[nx - 1][2][2] + N[nx - 1][2][3] + N[nx - 1][2][4] + N[nx - 1][2][6] + N[nx - 1][2][8]));

		N[nx - 1][ny - 3][3] = N[nx - 1][ny - 3][1]; N[nx - 1][ny - 3][4] = N[nx - 1][ny - 3][2]; N[nx - 1][ny - 3][7] = N[nx - 1][ny - 3][5];
		N[nx - 1][ny - 3][6] = N[nx - 1][ny - 3][8] = 0.5 * (rhoOutlet - (N[nx - 1][ny - 3][0] + N[nx - 1][ny - 3][1] + N[nx - 1][ny - 3][2] + N[nx - 1][ny - 3][3] + N[nx - 1][ny - 3][4] + N[nx - 1][ny - 3][5] + N[nx - 1][ny - 3][7]));

		//�����ܶȷֲ��ֽ�ֲ�����
		for (int j = 0; j < ny; j++) {
			if (Region::region[nx - 1][j] == Region::outlet || Region::region[nx - 1][j] == Region::FB_outlet)
				for (int k = 0; k < Q; k++) {
					N_r[nx - 1][j][k] = N[nx - 1][j][k] * rho_r[nx - 1][j] / (rho_r[nx - 1][j] + rho_b[nx - 1][j]);/////=================-------------===============
					N_b[nx - 1][j][k] = N[nx - 1][j][k] - N_r[nx - 1][j][k];
				}
		}

		return;
	}
}

double nu_(int& i, int j) {
	return (rho_r[i][j] + rho_b[i][j]) / (rho_r[i][j] / nu_r + rho_b[i][j] / nu_b);
}

double B(int& k) {
	if (k == 0) {
		return -4.0 / 27.0;
	}
	else if (k == 1 || k == 2 || k == 3 || k == 4) {
		return 2.0 / 27.0;
	}
	else {
		return 5.0 / 108.0;
	}
}

double omega(int& i, int& j) {
	return 2.0 / (6.0 * nu_(i, j) + 1.0);
}

double cosphi(int& k, double& Fc, double& FF) {
	if (FF > 0) {
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
	else
		return 0;
}

double cosphi_e(int& k, double& Fc, double& FF) {
	if (FF > 0.0) {
		//��һ�������������غ�������Ҫ,һ��Ҫ������������������������ɫ
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
	else
		return 0;
}

void streaming(double(&f)[nx][ny][Q]) {
	//������
	//����Ϊ streaming(N_r);
	for (int j = 0; j < ny; j++) {
		for (int i = nx - 1; i > 0; i--) {
			f[i][j][1] = f[i - 1][j][1];
		}
	}
	for (int i = 0; i < nx; i++) {
		for (int j = ny - 1; j > 0; j--) {
			f[i][j][2] = f[i][j - 1][2];
		}
	}
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx - 1; i++) {
			f[i][j][3] = f[i + 1][j][3];
		}
	}
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny - 1; j++) {
			f[i][j][4] = f[i][j + 1][4];
		}
	}
	for (int j = ny - 1; j > 0; j--) {
		for (int i = nx - 1; i > 0; i--) {
			f[i][j][5] = f[i - 1][j - 1][5];
		}
	}
	for (int j = ny - 1; j > 0; j--) {
		for (int i = 0; i < nx - 1; i++) {
			f[i][j][6] = f[i + 1][j - 1][6];
		}
	}
	for (int j = 0; j < ny - 1; j++) {
		for (int i = 0; i < nx - 1; i++) {
			f[i][j][7] = f[i + 1][j + 1][7];
		}
	}
	for (int j = 0; j < ny - 1; j++) {
		for (int i = nx - 1; i > 0; i--) {
			f[i][j][8] = f[i - 1][j + 1][8];
		}
	}
	return;
}

void out_file(const string& str,			//strΪ�ļ���
	const string& scalar_attribute_1,		//scalar_attribute_1Ϊ�������Ե���������
	const string& scalar_attribute_2,
	const string& scalar_attribute_3,
	const string& vector_attribute_1,		//vector_attribute_1Ϊʸ�����Ե���������
	const double(&u_vector_x)[nx][ny],		//u_vector_xΪ�ٶȳ���x����ķ���
	const double(&v_vector_y)[nx][ny],		//v_vector_yΪ�ٶȳ���y����ķ���
	const double(&Phi)[nx][ny],
	const areatype(&area)[nx][ny],
	const double(&rho)[nx][ny]) {
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
	outfile << "SCALARS" << " " << scalar_attribute_3 << " " << "double 1" << endl;
	outfile << "LOOKUP_TABLE default" << endl;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			outfile << rho[i][j] << endl;
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

void cal_Phi(double(&rho_r)[nx][ny], double(&rho_b)[nx][ny], double(&rho)[nx][ny]) {
	//�����ೡ
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] != Region::SL && Region::region[i][j] != Region::SB) {

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
	//2 4 6
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
			if (Region::region[i][j] != Region::SL && Region::region[i][j] != Region::SB) {

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
					//���������ص�omega(i,j)����ü�����ȶ�==============!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
	Region::setProcess();//����region
	Region::cal_ns();//�������߽�ĵ�λ��������ָ�����
	inoutBoundary::caculate_uInletDistribution();//������ڵ��ٶȷֲ�

	//��ʼ��
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] == Region::SL || Region::region[i][j] == Region::SB) {
				area[i][j] = solid;
			}
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] != Region::SL && Region::region[i][j] != Region::SB) {
				area[i][j] = blue;
			}
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (area[i][j] == blue && i < 5) {
				area[i][j] = red;
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
			if (Region::region[i][j] != Region::SL && Region::region[i][j] != Region::SB) {
				rho[i][j] = rho_r[i][j] + rho_b[i][j];
			}
		}
	}

	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] != Region::SL && Region::region[i][j] != Region::SB)
				for (int k = 0; k < Q; k++) {
					N_r[i][j][k] = w[k] * rho_r[i][j];//��ʼ��
					N_b[i][j][k] = w[k] * rho_b[i][j];
					N[i][j][k] = w[k] * rho[i][j];
				}
		}
	}

	double temp = 0;
	//��ʼ
	for (int step = 0; step <= 4000; step++) {
		if (step % 10 == 0) {
			temp = 0;
			cout << "��" << step << "��";
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					uu = u[i][j] * u[i][j] + v[i][j] * v[i][j];
					if (temp < sqrt(uu)) {
						temp = sqrt(uu);
					}
				}
			}
			printf("����ٶ� = %f\n", temp);
		}
		//������ײ����
		cal_Phi(rho_r, rho_b, rho);//������ɫ��,��׽����
		// �������ӿռ�ת�����ؿռ䡢�ؿռ����ײ���ؿռ�ת�������ӿռ�

		collision_single_phase(N, M, InvM, S, temp_moment, temp_moment_post_collison, rho, u, v);

		//����1��ʼ==================================================
		//debug(1, N_r, N_b, step);
		//����1����==================================================

		//�ೡ�ݶ�(�Ǳ߽���ೡ�ݶ�)
		for (int i = 0; i < nx; i++) {
			for (int j = 1; j < ny - 1; j++) {
				if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::inlet || Region::region[i][j] == Region::outlet) {

					if (i != 0 && i != nx - 1) {
						//������ɫ�ݶ�
						F_x[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i - 1][j + 1])
							+ 2.0 / 3.0 * (Phi[i + 1][j] - Phi[i - 1][j]) + 1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i - 1][j - 1]));

						F_y[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i + 1][j - 1])
							+ 2.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i - 1][j + 1] - Phi[i - 1][j - 1]));
					}
					else if (i == 0) {
						F_x[i][j] = 2.0 / 3.0 * (Phi[i + 1][j] - Phi[i][j]) + 1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i][j + 1]) +
							1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i][j - 1]);

						F_y[i][j] = 1.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i + 1][j - 1]);
					}
					else if (i == nx - 1) {
						F_x[i][j] = 2.0 / 3.0 * (Phi[i][j] - Phi[i - 1][j]) + 1.0 / 6.0 * (Phi[i][j + 1] - Phi[i - 1][j + 1]) +
							1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i][j - 1]);

						F_y[i][j] = 1.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i - 1][j + 1] - Phi[i - 1][j - 1]);
					}

					//��������1 ���ೡ�ݶȵĴ�С����aqrt(delta2)ʱ��Ϊ����
					if (F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j] > delta2) {
						area[i][j] = interface;
					}
					else if (Phi[i][j] > 0) {
						area[i][j] = red;
					}
					else
						area[i][j] = blue;
				}
			}
		}

		//ʪ��߽�����(�߽���ೡ�ݶ�)
		wettingBoundary::wetting_boundary(F_x, F_y);

		//����2��ʼ==================================================
		//debug(2, N_r, N_b, step);
		//����2����==================================================

		//�Ŷ�����
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (area[i][j] == interface || area[i][j] == contactline) {
					FF = F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j];
					for (int k = 0; k < Q; k++) {
						Fc = F_x[i][j] * cx[k] + F_y[i][j] * cy[k];
						N[i][j][k] = N[i][j][k] + A / 2.0 * sqrt(FF) * (w[k] * Fc * Fc / FF - B(k));
					}
				}
			}
		}

		//����3��ʼ==================================================
		//debug(3, N_r, N_b, step);
		//����3����==================================================

		//����ɫ����

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (Region::region[i][j] != Region::SL && Region::region[i][j] != Region::SB) {
					FF = F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j];
					for (int k = 0; k < Q; k++) {
						Fc = F_x[i][j] * cx[k] + F_y[i][j] * cy[k];

						N_r[i][j][k] = rho_r[i][j] / rho[i][j] * N[i][j][k] +
							beta * rho_r[i][j] * rho_b[i][j] / rho[i][j] * w[k] * cosphi_e(k, Fc, FF);
						N_b[i][j][k] = rho_b[i][j] / rho[i][j] * N[i][j][k] -
							beta * rho_r[i][j] * rho_b[i][j] / rho[i][j] * w[k] * cosphi_e(k, Fc, FF);

					}
				}
			}
		}

		//����4��ʼ==================================================
		//debug(4, N_r, N_b, step);
		//����4����==================================================

		//��
		streaming(N_r);
		streaming(N_b);

		//����5��ʼ==================================================
		//debug(5, N_r, N_b, step);
		//����5����==================================================

		//�߽�
		bounceBack(Region::region, N_r, u, v);
		bounceBack(Region::region, N_b, u, v);

		//ʹ���ֲܷ�������������ڳ��ڱ߽�

		inoutBoundary::inletVelocity();
		//inoutBoundary::inletPressure();
		//inoutBoundary::outletPressure();
		inoutBoundary::ouletNeumann(N_r);
		inoutBoundary::ouletNeumann(N_b);

		//����6��ʼ==================================================
		//debug(6, N_r, N_b, step);
		//����6����==================================================

		//�ܶ�
		for (int i = 0; i < nx; ++i) {
			for (int j = 0; j < ny; ++j) {
				if (Region::region[i][j] != Region::SL && Region::region[i][j] != Region::SB) {
					rho_r[i][j] = N_r[i][j][0] + N_r[i][j][1] + N_r[i][j][2] + N_r[i][j][3] +
						N_r[i][j][4] + N_r[i][j][5] + N_r[i][j][6] + N_r[i][j][7] + N_r[i][j][8];
					rho_b[i][j] = N_b[i][j][0] + N_b[i][j][1] + N_b[i][j][2] + N_b[i][j][3] +
						N_b[i][j][4] + N_b[i][j][5] + N_b[i][j][6] + N_b[i][j][7] + N_b[i][j][8];
					rho[i][j] = rho_r[i][j] + rho_b[i][j];
				}
			}
		}

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (Region::region[i][j] != Region::SL && Region::region[i][j] != Region::SB)
					for (int k = 0; k < Q; k++) {
						N[i][j][k] = N_r[i][j][k] + N_b[i][j][k];
					}
			}
		}

		//�ٶ�
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (Region::region[i][j] != Region::SL && Region::region[i][j] != Region::SB) {
					u[i][j] = (N[i][j][1] + N[i][j][5] + N[i][j][8] - N[i][j][3] - N[i][j][6] - N[i][j][7]) / rho[i][j];
					v[i][j] = (N[i][j][2] + N[i][j][5] + N[i][j][6] - N[i][j][4] - N[i][j][7] - N[i][j][8]) / rho[i][j];
				}
			}
		}

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (Region::region[i][j] == Region::SB || Region::region[i][j] == Region::SL) {
					Phi[i][j] = 0;
				}
			}
		}
		if (step % interval == 0) {
			out_file("new" + to_string(step / interval) + ".vtk", "Phi", "area", "density", "velocity", u, v, Phi, area, rho_r);
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
