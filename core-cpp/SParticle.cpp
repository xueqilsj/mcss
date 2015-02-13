/*
 * SParticle.cpp
 *
 *  Created on: May 20, 2009
 */
#include "SParticle.h"
#include <iostream> //< for cout
#include <cmath>
#include <float.h>
#ifndef M_PI
#define M_PI 3.1415926535897931
#endif

using namespace std;
SParticle::SParticle() {
	MP_XL = -2.0;
	MP_XH = 2.0;

	PAK_XL = -2.2;
	PAK_XH = 2.2;
}

SParticle::~SParticle() {
}

void SParticle::InitParameters(double initX, double off, double t) {
	uCurMove = 0;
	for (int i = 0; i < AI_NUM; i++)
		Arg_io[i] = 0;
	dX = initX;
	dOffset = off;
	Arg_io[AI_TEMPER] = t;
}
void SParticle::Gauge_MP() {
	// custom function F(x)
	if (dX < MP_XL || dX > MP_XH) {
		Arg_io[AI_HAMI] = DBL_MAX;
	} else if (-2 <= dX && dX <= -1.25) {
		Arg_io[AI_HAMI] = (1 + sin(2 * M_PI * dX));
	} else if (-1.25 <= dX && dX <= -0.25) {
		Arg_io[AI_HAMI] = 2 * (1 + sin(2 * M_PI * dX));
	} else if (-0.25 <= dX && dX <= 0.75) {
		Arg_io[AI_HAMI] = 3 * (1 + sin(2 * M_PI * dX));
	} else if (0.75 <= dX && dX <= 1.75) {
		Arg_io[AI_HAMI] = 4 * (1 + sin(2 * M_PI * dX));
	} else if (1.75 <= dX && dX <= 2) {
		Arg_io[AI_HAMI] = 5 * (1 + sin(2 * M_PI * dX));
	}
}

void SParticle::MetroplisTrial_MP(int Moves) {
	int CurStep = 0;
	double x1;
	double x2;

	while (CurStep < Moves) {
		//Delt();
		x1 = dX - dOffset;
		x2 = dX + dOffset;

		if ((MP_XL - x1) > 0) {
			x2 += (MP_XL - x1);
			x1 = MP_XL;
		}
		if ((x2 - MP_XH) > 0) {
			x1 -= (x2 - MP_XH);
			x2 = MP_XH;
		}

		dLastTrial = dX;
		dX = ran.Real(x1, x2);

		Arg_io[AI_DELT_HAMI] = Arg_io[AI_HAMI];//old
		Gauge_MP();//new
		Arg_io[AI_DELT_HAMI] = Arg_io[AI_HAMI] - Arg_io[AI_DELT_HAMI];
		if (Arg_io[AI_DELT_HAMI] > 0) {
			if (ran.Real() > exp(-1 * Arg_io[AI_DELT_HAMI] / Arg_io[AI_TEMPER])) {//random double value in the interval [0,1]
				//refused
				dX = dLastTrial;
				Arg_io[AI_HAMI] -= Arg_io[AI_DELT_HAMI];
			}

		} //else accepted

		uCurMove++;
		CurStep += 1;
	}
}

void SParticle::Gauge_Pak() {
	if (dX < PAK_XL || dX > PAK_XH) {
		Arg_io[AI_HAMI] = DBL_MAX;
	} else {
		Arg_io[AI_HAMI] = ((dX + 1) * (dX + 1) - 1) * ((dX - 1) * (dX - 1)
				- 0.9);
	}
}
void SParticle::MetroplisTrial_Pak(int Moves) {

	int CurStep = 0;
	double x1;
	double x2;
	while (CurStep < Moves) {
		//Delt();
		x1 = dX - dOffset;
		x2 = dX + dOffset;

		if ((PAK_XL - x1) > 0) {
			x2 += (PAK_XL - x1);
			x1 = PAK_XL;
		}
		if ((x2 - PAK_XH) > 0) {
			x1 -= (x2 - PAK_XH);
			x2 = PAK_XH;
		}

		dLastTrial = dX;
		dX = ran.Real(x1, x2);

		Arg_io[AI_DELT_HAMI] = Arg_io[AI_HAMI];//old
		Gauge_Pak();//new
		Arg_io[AI_DELT_HAMI] = Arg_io[AI_HAMI] - Arg_io[AI_DELT_HAMI];

		if (Arg_io[AI_DELT_HAMI] > 0) {
			if (ran.Real() > exp(-1 * Arg_io[AI_DELT_HAMI] / Arg_io[AI_TEMPER])) {//random double value in the interval [0,1]
				//refused
				dX = dLastTrial;
				Arg_io[AI_HAMI] -= Arg_io[AI_DELT_HAMI];
			}

		} //else accepted

		uCurMove++;
		CurStep += 1;
	}
}
void SParticle::SetArg(int site, double df) {
	if (site >= 0 && site < AI_NUM)
		Arg_io[site] = df;
}
double SParticle::GetArg(int site) {
	if (site >= 0 && site < AI_NUM)
		return Arg_io[site];
	else {
		cout << "#Boundary error:" << site << endl;
		return 0;
	}
}
