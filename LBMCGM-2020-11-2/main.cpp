#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "MRT.h"
//根据文献"Lattice Boltzmann simulation of pressure-driven two-phase flows in capillary tube and porous medium"
//以及huang2014的文章对出入口边界进行调整
//对于侵入的流体(红色):入口为速度边界，出口为压力为0的边界
//对于抗侵入的流体(蓝色):入口为压力为0的边界，出口为恒压的边界
//具体而言：   
//		红色 :	入口 : 速度为 u0
//				出口 : 压力为 极小值 例如
//		蓝色 :  入口与出口均为压力极小值
//修改了sourceTermPrecondition()函数，将不在界面的表面张力力赋值为0

//根据文献"Lattice Boltzmann simulation of pressure-driven two-phase flows in capillary tube and porous medium"
//从总分布函数的角度设计入口出口边界

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
double rhoIn = 1e-3;//这是为蓝色流体准备的入口0压力近似值
double rhoOut = 0.1;//对于在出口处的流体的密度
double u[nx][ny]{ 0 };
double v[nx][ny]{ 0 };
double uInlet = 1e-3;//入口流速
double uInletDistribution[ny];
double vInlet = 0.0;
//注意，调整uInlet rhoOutl的大小以保证最大流速在0.1以下。（前1000步必须在0.1以下才能保证后续的稳定）
double Phi[nx][ny]{ 0 };//相场
double F_x[nx][ny]{ 0 };//相场梯度
double F_y[nx][ny]{ 0 };
double n[nx][ny][2]{ 0 };//相场梯度的单位向量
double kappa[nx][ny]{ 0 };//曲率
double F[nx][ny][2]{ 0 };//外力，这里用来描述表面张力力，也可以加入其它体力
double cx[9]{ 0,1,0,-1,0,1,-1,-1,1 };
double cy[9]{ 0,0,1,0,-1,1,1,-1,-1 };
double w[9] = { 4.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 9.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0,1.0 / 36.0 };
double rho_r0 = 1;//初始密度
double rho_b0 = 1;
double nu_r = 1.0 / 6.0;//粘滞系数，粘性越大，弛豫率就越大,趋于平衡的速度也就越快,这个也会影响伪电流，两个流体的粘性都变大时，伪电流会减小
double nu_b = 1.0 / 60.0;
double contactAngle = 120.0 / 180.0 * 3.141592653;//接触角!!!!!!!!!!!!!!!!!
double theta = 0.004;//表面张力,这个参数如果较大时，边界会不稳定,这是个很重要的参数，界面受力由它控制，伪电流也受其影响，经过测试，最好不大于0.005
double beta = 0.99;
double delta2 = 0.05;//根据相场梯度控制界面厚度,目前这个值得取法是根据图像试出来的
double uu, cu, FF, Fc;
// S对角矩阵 松弛系数矩阵 s7与s8与BGK中的omega相同 
double s0 = 1.0, s1 = 1.63, s2 = 1.14, s3 = 1.0, s4 = 1.92, s5 = 1.0, s6 = 1.92, s7, s8;//s7 s8会变化
double S[9] = { s0,s1,s2,s3,s4,s5,s6,s7,s8 };
// 平衡矩分布数组
double temp_moment[9]{ 0 };
double temp_moment_eq[9]{ 0 };
double temp_moment_post_collison[9]{ 0 };
double temp_moment_F[9]{ 0 };
double temp_moment_post_collison_F[9]{ 0 };
enum areatype { interface, red, blue, other };
areatype area[nx][ny];

int interval = 10;

namespace Region {
	enum regiontype { F, FB, FL, S, SB, SL };//一般是FB,FL,SB,SL
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
		//修正ns
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
		//用setS()函数设计好固体区域后，直接调用Region::setRegion::setProcess();
		void setS() {
			//只要设计好这个函数就可以了

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
		//湿润边界条件,这个边界条件只使用于两相界面
		//在边界使用虚拟的梯度，只在同时在接触线和边界上的晶格设置
		double temp1, temp2;
		//第一步先估算SB或者S的相场Phi
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
		//第二步利用第一步估算的SB处的相场Phi计算FB处的相场梯度F_x,F_y
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				if (Region::region[i][j] == Region::FB) {

					F_x[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i - 1][j + 1])
						+ 2.0 / 3.0 * (Phi[i + 1][j] - Phi[i - 1][j]) + 1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i - 1][j - 1]));

					F_y[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i + 1][j - 1])
						+ 2.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i - 1][j + 1] - Phi[i - 1][j - 1]));

					//划分区域2 当相场梯度的大小大于aqrt(delta2)时认为界面
					if (F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j] > delta2) {
						area[i][j] = interface;
					}
					else {
						area[i][j] = other;
					}
				}
			}
		}
		//第三步保持在FB估算的相场梯度F_x,F_y大小不变，根据接触角与FB处的固壁法向ns修正方向
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

		//合并分布函数
		for (int j = 2; j <= ny - 3; j++) {
			for (int k = 0; k < Q; k++) {
				f[0][j][k] = f_r[0][j][k] + f_b[0][j][k];
			}
		}

		//入口的单位法向为(1,0)时的速度边界，(1,0)指向里面
		for (int j = 3; j <= ny - 4; j++) {
			rho[0][j] = ((f[0][j][0] + f[0][j][2] + f[0][j][4]) + 2.0 * (f[0][j][3] + f[0][j][6] + f[0][j][7])) / (1 - u[j]);
			f[0][j][1] = f[0][j][3] + 2.0 * rho[0][j] * u[j] / 3.0;
			f[0][j][5] = f[0][j][7] - 0.5 * (f[0][j][2] - f[0][j][4]) + rho[0][j] * u[j] / 6.0 + 0.5 * rho[0][j] * vInlet;
			f[0][j][8] = f[0][j][6] + 0.5 * (f[0][j][2] - f[0][j][4]) + rho[0][j] * u[j] / 6.0 - 0.5 * rho[0][j] * vInlet;
		}

		//入口的corner处需要特别处理
		//左下角
		f[0][2][1] = f[0][2][3]; f[0][2][2] = f[0][2][4]; f[0][2][5] = f[0][2][7];
		f[0][2][6] = f[0][2][8] = 0.5 * (rho[0][3] - (f[0][2][0] + f[0][2][1] + f[0][2][2] + f[0][2][3] + f[0][2][4] + f[0][2][5] + f[0][2][7]));
		//左上角
		f[0][ny - 3][1] = f[0][ny - 3][3]; f[0][ny - 3][4] = f[0][ny - 3][2]; f[0][ny - 3][8] = f[0][ny - 3][6];
		f[0][ny - 3][5] = f[0][ny - 3][7] = 0.5 * (rho[0][ny - 4] - (f[0][ny - 3][0] + f[0][ny - 3][1] + f[0][ny - 3][2] + f[0][ny - 3][3] + f[0][ny - 3][4] + f[0][ny - 3][6] + f[0][ny - 3][8]));

		//根据密度分布分解分布函数
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
		//这个函数分别调用不同的入口密度
		//入口的单位法向为(1,0)时的压力边界,(1,0)指向里面

		//合并分布函数
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

		//出口的corner处需要特别处理
		//左上角
		f[0][ny - 3][1] = f[0][ny - 3][3];
		f[0][ny - 3][4] = f[0][ny - 3][2];
		f[0][ny - 3][8] = f[0][ny - 3][6];
		f[0][ny - 3][5] = f[0][ny - 3][7] = 0.5 * (rhoIn - (f[0][ny - 3][1] + f[0][ny - 3][2] +
			f[0][ny - 3][3] + f[0][ny - 3][4] + f[0][ny - 3][6] + f[0][ny - 3][8]));
		//左下角
		f[0][2][1] = f[0][2][3];
		f[0][2][5] = f[0][2][7];
		f[0][2][4] = f[0][2][2];
		f[0][2][6] = f[0][2][8] = 0.5 * (rhoIn - (f[0][2][1] + f[0][2][2] +
			f[0][2][3] + f[0][2][4] + f[0][2][5] + f[0][2][7]));

		//根据密度分布来调整这里的分布函数
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
		//这个函数分别调用不同的出口密度
		//出口的单位法向为(1,0)时的压力边界,(1,0)指向外面

		//合并分布函数
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

		//出口的corner处需要特别处理
		//右上角
		f[nx - 1][ny - 3][3] = f[nx - 1][ny - 3][1];
		f[nx - 1][ny - 3][7] = f[nx - 1][ny - 3][5];
		f[nx - 1][ny - 3][4] = f[nx - 1][ny - 3][2];
		f[nx - 1][ny - 3][6] = f[nx - 1][ny - 3][8] = 0.5 * (rhoOut - (f[nx - 1][ny - 3][1] + f[nx - 1][ny - 3][2] +
			f[nx - 1][ny - 3][3] + f[nx - 1][ny - 3][4] + f[nx - 1][ny - 3][5] + f[nx - 1][ny - 3][7]));
		//右下角
		f[nx - 1][2][2] = f[nx - 1][2][4];
		f[nx - 1][2][3] = f[nx - 1][2][1];
		f[nx - 1][2][6] = f[nx - 1][2][8];
		f[nx - 1][2][5] = f[nx - 1][2][7] = 0.5 * (rhoOut - (f[nx - 1][2][1] + f[nx - 1][2][2] +
			f[nx - 1][2][3] + f[nx - 1][2][4] + f[nx - 1][2][6] + f[nx - 1][2][8]));

		//根据密度分布来调整这里的分布函数
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
	//流函数
	//调用为 streaming(N_r);
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

void out_file(const string& str,			//str为文件名
	const string& scalar_attribute_1,		//scalar_attribute_1为标量属性的量的名称
	const string& vector_attribute_1,		//vector_attribute_1为矢量属性的量的名称
	const double(&u_vector_x)[nx][ny],		//u_vector_x为速度场的x方向的分量
	const double(&v_vector_y)[nx][ny],		//v_vector_y为速度场的y方向的分量
	const double(&Phi)[nx][ny]) {			//rho_scalar_1为相场
	//后处理文件输出函数
	//这是多相流的后处理输出的文件格式
	//调用为：
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

void out_file(const string& str,			//str为文件名
	const string& scalar_attribute_1,		//scalar_attribute_1为标量属性的量的名称
	const string& scalar_attribute_2,
	const string& vector_attribute_1,		//vector_attribute_1为矢量属性的量的名称
	const double(&u_vector_x)[nx][ny],		//u_vector_x为速度场的x方向的分量
	const double(&v_vector_y)[nx][ny],		//v_vector_y为速度场的y方向的分量
	const double(&Phi)[nx][ny],
	const areatype(&area)[nx][ny]) {
	//后处理文件输出函数
	//这是多相流的后处理输出的文件格式
	//调用为：
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
	//调用为 bounceBack(Region::region,f,u,v);
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
	//相场界面的分散化
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
	//计算颜色场，并划分区域
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
	//计算临时平衡矩
	//形参rho对应实参rho[i][j]
	//形参u对应实参u[i][j]
	//形参v对应实参v[i][j]
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
	// 将从粒子空间转换到矩空间、矩空间的碰撞、矩空间转换回粒子空间
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {

				// 转换到矩空间
				for (int k = 0; k < 9; k++) {
					temp_moment[k] = M[k][0] * f[i][j][0] + M[k][1] * f[i][j][1] + M[k][2] * f[i][j][2] +
						M[k][3] * f[i][j][3] + M[k][4] * f[i][j][4] + M[k][5] * f[i][j][5] +
						M[k][6] * f[i][j][6] + M[k][7] * f[i][j][7] + M[k][8] * f[i][j][8];
				}

				//计算临时平衡矩
				calculate_temp_moment_eq(rho[i][j], u[i][j], v[i][j]);//这里的rho为函数的参数而不是全局变量

				// 矩空间的碰撞
				//S[7]与S[8]在界面处是变化的
				for (int k = 0; k < 7; k++) {
					temp_moment_post_collison[k] = (1 - S[k]) * temp_moment[k] + S[k] * temp_moment_eq[k];
					//temp_moment_post_collison[k] = (1 - 0.8 * omega(i, j)) * temp_moment[k] + 0.8 * omega(i, j) * temp_moment_eq[k];
				}
				temp_moment_post_collison[7] = (1 - omega(i, j)) * temp_moment[7] + omega(i, j) * temp_moment_eq[7];
				temp_moment_post_collison[8] = (1 - omega(i, j)) * temp_moment[8] + omega(i, j) * temp_moment_eq[8];

				// 转换回粒子空间,分布函数 f 经过了矩空间内的碰撞
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
	//形参F[2]对应实参F[i][j]
	//形参u对应实参u[i][j]
	//形参v对应实参v[i][j]
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
	//使用该源项给流体施加表面张力 (预处理)

	//第一步计算相场梯度单位向量
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

	//第二步计算 曲率 第三步 表面张力力 第四步 源项 MF
	for (int i = 1; i < nx - 1; i++) {
		for (int j = 1; j < ny - 1; j++) {
			if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
				if (area[i][j] == interface) {
					//曲率
					kappa[i][j] = -(1.0 / 12.0 * (n[i + 1][j + 1][0] - n[i - 1][j + 1][0] + n[i + 1][j - 1][0] - n[i - 1][j - 1][0] +
						n[i + 1][j + 1][1] - n[i + 1][j - 1][1] + n[i - 1][j + 1][1] - n[i - 1][j - 1][1]) +
						1.0 / 3.0 * (n[i + 1][j][0] - n[i - 1][j][0] + n[i][j + 1][1] - n[i][j - 1][1]));

					//表面张力力
					F[i][j][0] = 0.5 * theta * kappa[i][j] * F_x[i][j];
					F[i][j][1] = 0.5 * theta * kappa[i][j] * F_y[i][j];
					//源项 M^(-1) * (I  0.5 * S) * MF 
					//放在main中计算
				}
				else
				{
					//注意将不在界面的表面张力力置为0
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
	Region::setRegion::setProcess();//设置region
	Region::cal_ns();//计算固体边界的单位法向量，指向固体
	inoutBoundary::caculate_uInletDistribution();//计算入口的速度分布
	//初始化
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

	//设计液滴的位置参数
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
				/*ppt中使用的是这个初始分布!!!!!!!!!!!!!!!!!
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
				N_r[i][j][k] = w[k] * rho_r[i][j];//初始化
				N_b[i][j][k] = w[k] * rho_b[i][j];
			}
		}
	}
	double temp = 0;
	//开始
	for (int step = 0; step < 8001; step++) {
		if (step % 10 == 0) {
			temp = 0;
			//cout << "第" << step << "步" << endl;
			cout << "第" << step << "步";
			for (int i = 0; i < nx; i++) {
				for (int j = 0; j < ny; j++) {
					uu = u[i][j] * u[i][j] + v[i][j] * v[i][j];
					if (temp < sqrt(uu)) {
						temp = sqrt(uu);
					}
				}
			}
			cout << "最大速度为:" << temp << endl;
		}
		//单相碰撞算子
		cal_Phi(rho_r, rho_b, rho);//计算颜色场,捕捉界面
		// 将从粒子空间转换到矩空间、矩空间的碰撞、矩空间转换回粒子空间

		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				for (int k = 0; k < Q; k++) {
					N[i][j][k] = N_r[i][j][k] + N_b[i][j][k];
				}
			}
		}

		collision_single_phase(N, M, InvM, S, temp_moment, temp_moment_post_collison, rho, u, v);

		//调试1开始==================================================
		debug(1, N_r, N_b, step);
		//调试1结束==================================================

		//扰动算子
		for (int i = 1; i < nx - 1; i++) {
			for (int j = 1; j < ny - 1; j++) {
				if (Region::region[i][j] == Region::FL) {

					//计算颜色梯度
					F_x[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i - 1][j + 1])
						+ 2.0 / 3.0 * (Phi[i + 1][j] - Phi[i - 1][j]) + 1.0 / 6.0 * (Phi[i + 1][j - 1] - Phi[i - 1][j - 1]));
					F_y[i][j] = 0.5 * (1.0 / 6.0 * (Phi[i + 1][j + 1] - Phi[i + 1][j - 1])
						+ 2.0 / 3.0 * (Phi[i][j + 1] - Phi[i][j - 1]) + 1.0 / 6.0 * (Phi[i - 1][j + 1] - Phi[i - 1][j - 1]));

					//划分区域1 当相场梯度的大小大于aqrt(delta2)时认为界面
					if (F_x[i][j] * F_x[i][j] + F_y[i][j] * F_y[i][j] > delta2) {
						area[i][j] = interface;
					}
					else {
						area[i][j] = other;
					}
				}
			}
		}

		//湿润边界条件

		wettingBoundary::wetting_boundary(F_x, F_y);

		//调试2开始==================================================
		debug(2, N_r, N_b, step);
		//调试2结束==================================================

		//计算表面张力力（源项的预处理）
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

						// 转换回粒子空间,分布函数 f 经过了矩空间内的碰撞
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

		//调试3开始==================================================
		debug(3, N_r, N_b, step);
		//调试3结束==================================================

		//重着色算子

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (Region::region[i][j] == Region::FL || Region::region[i][j] == Region::FB) {
					if (area[i][j] == interface) {
						//一般来说，这里应该添加一个判断，仅在interface处使用
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

		//调试4开始==================================================
		debug(4, N_r, N_b, step);
		//调试4结束==================================================

		//流
		streaming(N_r);
		streaming(N_b);

		//调试5开始==================================================
		debug(5, N_r, N_b, step);
		//调试5结束==================================================

		//边界
		bounceBack(Region::region, N_r, u, v);
		bounceBack(Region::region, N_b, u, v);

		//使用总分布函数来计算入口出口边界
		inoutBoundary::inletVelocity(N, N_r, N_b, rho, rho_r, rho_b, u, v, uInletDistribution);
		inoutBoundary::outletPressure(N, N_r, N_b, rho, rho_r, rho_b, u, v, rhoOut);



		//调试6开始==================================================
		debug(6, N_r, N_b, step);
		//调试6结束==================================================

		//密度
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

		//速度
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
		//调试7开始==================================================
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
		//调试7结束==================================================
	}
	return 0;
}