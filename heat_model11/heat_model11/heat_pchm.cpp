#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include "func_header.h"
using namespace std;

const int Nr = 100;					//	size in r-direction
const int Nz = 54;					//	size in z-direction
double T[Nr][Nz], T_next[Nr][Nz];	//	temperature in current and next moment
double F[Nr][Nz], F_next[Nr][Nz], I[Nr][Nz];	//	reaction coefficient in current and next moment, inetnsity of laser beam
int sec = 0;						//  time (sec)
//output files
ofstream trz0("trz0.csv"), trzh("trzh.csv"), tz("tz.csv"), frz0("frz0.csv"), frzh("frzh.csv"), fz("fz.csv");

void heatingPCHM3D() {
				
	int time_total = 100;	//	total time in sec

	double lambda = 0.16;		//	thermal conductivity
	double C = 1240;			//	heat capacity
	double Lx = 5e-3;			//	glass radius
	double Lz = 2.7e-3;			//	depth
	double I0 = 1592;			//	peak value of laser's intensity
	double p = 1180;			//	density
	double b = 1.8e-5;			//	reaction's speed
	double R = 0.002;			//	laser's beam radius

	double a = lambda / (C*p);	//	coefficient a
	double h = 5.0e-5;			//	spatial step
	double k0 = 528, k1 = 214;	//	exponential transmittance coefficient of glass
	double alpha = 1 / pow(R, 2);	//	alpha coefficient
	double t = 0.001*pow(h, 2) / a;	//	time step

	//initializing initial values of intensity, temperature and concetration
	for (int i = 0; i < Nr; i++) {
		double r = i*h;
		I[i][0] = I0*exp(-alpha* pow(r, 2));
		for (int j = 0; j < Nz; j++) {
			T[i][j] = 0;
			F[i][j] = 1;
		}
	}

	//main calculation loop
	for (int count = 0; t*count <= time_total; count++) {

		//calculation laser intensity throughout glass
		for (int j = 0; j < Nz - 1; j++) {
			for (int i = 0; i < Nr; i++) {
				I[i][j + 1] = (-(k1 + (k0 - k1)*F[i][j]))*h*I[i][j] + I[i][j];
			}
		}
		//calculation temperature
		for (int j = 1; j < Nz - 1; j++) {
			for (int i = 1; i < Nr - 1; i++) {
				double r = i*h;
				double Q;
				Q = -(I[i][j + 1] - I[i][j]) / (h);
				T_next[i][j] = T[i][j] + a*t*((T[i + 1][j] + T[i - 1][j] + T[i][j + 1] + T[i][j - 1] - 4 * T[i][j]) / pow(h, 2.0) + (1 / r) * (T[i + 1][j] - T[i - 1][j]) / (2 * h)) + Q*t / (C*p);
				//boundary conditions
				if (j == 1) T_next[i][0] = lambda*T[i][1] / (lambda - 10.0*h);
				if (j == Nz - 2) T_next[i][Nz - 1] = lambda*T[i][Nz - 2] / (lambda + 10.0*h);
			}
			//boundary conditions
			T_next[0][j] = lambda*T_next[1][j] / (lambda - 10.0*h);
			T_next[Nr - 1][j] = lambda*T_next[Nr - 2][j] / (lambda + 10.0*h);
		}
		T_next[0][0] = T_next[1][0];

		for (int j = 0; j < Nz; j++) {
			for (int i = 0; i < Nr; i++) {
				F_next[i][j] = -b*t*I[i][j] * F[i][j] + F[i][j];	//reaction coefficient calculations
				T[i][j] = T_next[i][j];
				F[i][j] = F_next[i][j];
			}
		}

		if (count % 1000 == 0) cout << count << endl;	//just to check progress

		if (count*t >= sec) {
			sendData();
			++sec;
		}

	}
	tz.close(); fz.close(); trz0.close(); trzh.close(); frz0.close(); frzh.close();

}


void sendData() {
	tz << "t= " << sec << " ;";
	for (int i = 0; i < Nr; ++i) {
		trz0 << T[i][0] << ";";
		trzh << T[i][Nz - 1] << ";";
		frz0 << F[i][0] << ";";
		frzh << F[i][Nz - 1] << ";";
	}
	trz0 << endl;
	trzh << endl;
	frz0 << endl;
	frzh << endl;
	for (int j = 0; j < Nz; ++j) {
		tz << T[0][j] << ";";
		fz << F[0][j] << ";";

	}
	tz << endl;
	fz << endl;
}

