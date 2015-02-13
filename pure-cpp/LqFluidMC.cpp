#include "LqFluid.h"
#include <cmath>
using namespace std;
#include <iostream>
#include <iomanip>
#include <time.h>

#define DBL_EPSILON    2e-8
#define SWITCHDIST 8.5 
#define CUTOFF     10 
#define PAIRLISTDIST 11.5 

int main(int argc, char *argv[]) {

	int nsteps = 1000;

	int nl0 = 20, nl1 = 20;
	int natoms = nl0 * nl1;

	LqFluid lqfluid;

	lqfluid.InitParameters(0.87, 0.2);
	lqfluid.ResetLQ(119, 0.341, SWITCHDIST, CUTOFF);
	lqfluid.dBeta = 2.42;

	double *Corxs = new double[natoms];
	double *Corys = new double[natoms];
	double *Fxs = new double[natoms];
	double *Fys = new double[natoms];
	lqfluid.SetConf(Corxs, Corys, Fxs, Fys, nl0, nl1);
	lqfluid.InitConf(1);
	lqfluid.Gauge();
	cout << lqfluid.Evdw << endl;
	lqfluid.Arg_io[lqfluid.AI_HAMI] = lqfluid.Evdw;
	clock_t __start = clock();
	for (int i = 0; i < nsteps; i++) {
		lqfluid.MetroplisSweepBeta(10);
		//lqfluid.Gauge();
		cout << setprecision(16) << lqfluid.Arg_io[lqfluid.AI_DELT_HAMI]
				<< '\t' << lqfluid.Arg_io[lqfluid.AI_HAMI] << '\t'
				<< lqfluid.dMax << '\t' << lqfluid.dAcceptRatio << endl;
	}
	cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\tseconds";

	return 0;
}
