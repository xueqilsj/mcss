/*
 * Linear.cpp
 *
 *  Created on: May 20, 2009
 */
#include "Linear.h"
#include <iostream> //< for cout
#include <cmath>
#include <float.h>
#ifndef M_PI
#define M_PI 3.1415926535897931
#endif

using namespace std;
Linear::Linear() {
}

Linear::~Linear() {
}

void Linear::InitParameters(double initp, double off, double T) {
	uCurMove = 0;
	for (int i = 0; i < AI_NUM; i++)
		Arg_io[i] = 0;
	//TODO: depend on custom function F(x)
	dDim[0] = -2.0;
	dDim[1] = 2.0;
	dPoint[0] = initp;
	dOffset[0] = off;
	Arg_io[AI_TEMPER] = T;
}
void Linear::Gauge() {
	Arg_io[AI_MAG] = 0;
	//TODO: depend on custom function F(x)
	if (dPoint[0] < -2 || dPoint[0] > 2) {
		Arg_io[AI_HAMI] = DBL_MAX;
	} else if (-2 <= dPoint[0] && dPoint[0] <= -1.25) {
		Arg_io[AI_HAMI] = (1 + sin(2 * M_PI * dPoint[0]));
	} else if (-1.25 <= dPoint[0] && dPoint[0] <= -0.25) {
		Arg_io[AI_HAMI] = 2 * (1 + sin(2 * M_PI * dPoint[0]));
	} else if (-0.25 <= dPoint[0] && dPoint[0] <= 0.75) {
		Arg_io[AI_HAMI] = 3 * (1 + sin(2 * M_PI * dPoint[0]));
	} else if (0.75 <= dPoint[0] && dPoint[0] <= 1.75) {
		Arg_io[AI_HAMI] = 4 * (1 + sin(2 * M_PI * dPoint[0]));
	} else if (1.75 <= dPoint[0] && dPoint[0] <= 2) {
		Arg_io[AI_HAMI] = 5 * (1 + sin(2 * M_PI * dPoint[0]));
	}
}

double Linear::Delt() {
	double x1 = dPoint[0] - dOffset[0];
	double x2 = dPoint[0] + dOffset[0];

	if ((dDim[0] - x1) > 0) {
		x2 += (dDim[0] - x1);
		x1 = dDim[0];
	}
	if ((x2 - dDim[1]) > 0) {
		x1 -= (x2 - dDim[1]);
		x2 = dDim[1];
	}
	dLastTrial[0] = dPoint[0];
	dPoint[0] = ran.Real(x1, x2);

	Arg_io[AI_DELT_HAMI] = Arg_io[AI_HAMI];//old
	Gauge();//new
	Arg_io[AI_DELT_HAMI] = Arg_io[AI_HAMI]
			- Arg_io[AI_DELT_HAMI];
	return dLastTrial[0];
}

int Linear::Accepted() {//1 absolutely accepted,0 modestly accepted,-1  refused
	if (Arg_io[AI_DELT_HAMI] > 0) {
		Arg_io[AI_EXP_DELT_HAMI] = exp(-1 * Arg_io[AI_DELT_HAMI]
				/ Arg_io[AI_TEMPER]);
		Arg_io[AI_RANDOM] = ran.Real();//< random double value in the interval [0,1]
		if (Arg_io[AI_RANDOM] < Arg_io[AI_EXP_DELT_HAMI]) {
			nAccepted = 0;
		} else {//< refused
			Arg_io[AI_HAMI] -= Arg_io[AI_DELT_HAMI];
			//Arg_io[AI_MAG] -= Arg_io[AI_DELT_MAG];
			dPoint[0] = dLastTrial[0];
			nAccepted = -1;
		}
	} else {

		nAccepted = 1;
	}
	return nAccepted;
}

void Linear::MetroplisTrial(int Moves) {
	int CurStep = 0;
	while (CurStep < Moves) {
		Delt();
		Accepted();
		uCurMove++;
		CurStep += 1;
	}
}
void Linear::SetArg(int site, double df) {
	if (site >= 0 && site < AI_NUM)
		Arg_io[site] = df;
}
double Linear::GetArg(int site) {
	if (site >= 0 && site < AI_NUM)
		return Arg_io[site];
	else {
		cout << "#Boundary error:" << site << endl;
		return 0;
	}
}
double Linear::GetPoint(int in) {
	return dPoint[in];
}
