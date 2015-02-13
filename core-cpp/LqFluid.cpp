/*
 * LqFluid.cpp
 *  Basic C++ Class for LqFluid model
 *  Created on: Dec. 29, 2009
 */
#include "LqFluid.h"
#include <iomanip>
#include <iostream> //< for cout
#include <cmath> //exp,sqrt
#include <ctype.h> //isdigit
#define DBL_EPSILON    2e-8
using namespace std;

// Return the next nearest integer resulting from
// the divison of x by pbcx.(abs(x))

static int nint(double x, double pbcx) {
	int nInt;
	double div;
	div = x / pbcx;
	if (x >= 0)
		nInt = (int) (div + 0.5);
	else
		nInt = -(int) (0.5 - div);
	return nInt;
}

LqFluid::LqFluid() {
}

LqFluid::~LqFluid() {
}

void LqFluid::InitParameters(double ds, double dmax) {
	uCurMove = 0;
	uAccepted = 0;
	for (int i = 0; i < AI_NUM; i++)
		Arg_io[i] = 0;
	Arg_io[AI_TEMPER] = 1;
	dBeta = 1;
	Evdw = 0;

	dDist = ds;
	dMax = dmax;
}

void LqFluid::SetConf(double *Corxs, double *Corys, double *Fxs, double *Fys,
		int dim0, int dim1) {
	nDim0 = dim0;
	nDim1 = dim1;
	nN = nDim0 * nDim1;
	Xp = Corxs;
	Yp = Corys;
	FXp = Fxs;
	FYp = Fys;
	pbcx = dDist * nDim0;//?remain space
	pbcy = dDist * nDim1;
}

void LqFluid::InitConf(int d) {
	if (d > -2) {
		int nl0 = nDim0, nl1 = nDim1;
		for (int j = 0; j < nl1; j++) {
			for (int i = 0; i < nl0; i++) {
				Xp[nl0 * j + i] = dDist * (i - nl0 / 2);
				Yp[nl0 * j + i] = dDist * (j - nl1 / 2);
				//cout << Xp[nl1 * j + i] << ':' << Yp[nl1 * j + i] << endl;
				FXp[nl0 * j + i] = 0;
				FYp[nl0 * j + i] = 0;

			}

		}
	}
}
void LqFluid::ResetLQ(double epsil, double sig, double swit, double cutf) {
	cutoff = cutf;
	cut2 = cutoff * cutoff;
	switchdist = swit;
	switch2 = switchdist * switchdist;
	epsilon = epsil;
	sigma = sig;

	sig *= sig * sig;
	sig *= sig;//s^6
	vdwB = 4.0 * sig * epsil;
	vdwA = vdwB * sig;
}
double LqFluid::GetLqValue(double dist) {
	double AmBterm, r, r_1, r_2, r_6, r_12, switchVal, dSwitchVal, c1, c3, c2,
			c4;

	c1 = 1.0 / (cut2 - switch2);
	c1 = c1 * c1 * c1;
	c3 = 4 * c1;

	if (dist > cut2) {
		return 0;
	} else {
		r = sqrt(dist);
		r_1 = 1.0 / r;
		r_2 = r_1 * r_1;
		r_6 = r_2 * r_2 * r_2;
		r_12 = r_6 * r_6;
		switchVal = 1, dSwitchVal = 0;
		if (dist > switch2) {
			c2 = cut2 - dist;
			c4 = c2 * (cut2 + 2 * dist - 3.0 * switch2);
			switchVal = c2 * c4 * c1;
			dSwitchVal = c3 * r * (c2 * c2 - c4);
		}
		AmBterm = (vdwA * r_6 - vdwB) * r_6;
		return switchVal * AmBterm;
	}
}
int LqFluid::GetRDF(double *rdf, int len) {
	int ni = (nN - 1) * nN / 2;
	if (len != ni)
		return -1;

	ni = 0;
	for (unsigned int i = 0; i < nN - 1; i++) {
		for (unsigned int j = i + 1; j < nN; j++) {
			rdf[ni++] = sqrt(pbc_length_ji2(j, i));
		}
	}
	return ni;
}

double LqFluid::pbc_length_ji2(int j, int i) {
	//
	double rx = Xp[j] - Xp[i];
	double ry = Yp[j] - Yp[i];
	rx = Xp[i] + rx - pbcx * nint(rx, pbcx);// rint(rx / pbcx);
	ry = Yp[i] + ry - pbcy * nint(ry, pbcy);
	//cout << Xp[j] << '\t' << rx << "\tv" << endl;
	//cout << Yp[j] << '\t' << ry << "\t^" << endl;
	return ((Xp[i] - rx) * (Xp[i] - rx) + (Yp[i] - ry) * (Yp[i] - ry));
}
void LqFluid::addF_ji(int j, int i, double f) {
	double rx = Xp[j] - Xp[i];
	double ry = Yp[j] - Yp[i];
	rx = Xp[i] + rx - pbcx * nint(rx, pbcx);
	ry = Yp[i] + ry - pbcy * nint(ry, pbcy);

	FXp[i] -= (rx - Xp[i]) * f;
	FYp[i] -= (ry - Yp[i]) * f;

	FXp[j] += (rx - Xp[i]) * f;
	FYp[j] += (ry - Yp[i]) * f;
}
void LqFluid::SetArg(int site, double df) {
	if (site >= 0 && site < AI_NUM)
		Arg_io[site] = df;
}
double LqFluid::GetArg(int site) {
	if (site >= 0 && site < AI_NUM)
		return Arg_io[site];
	else {
		cout << "#Boundary error:" << site << endl;
		return 0;
	}
}

void LqFluid::Gauge() {
	double force_r, AmBterm, r, r_1, r_2, r_6, r_12, switchVal, dSwitchVal, c1,
			c3, c2, c4, dist;

	Evdw = 0;
	c1 = 1.0 / (cut2 - switch2);
	c1 = c1 * c1 * c1;
	c3 = 4 * c1;
	for (unsigned int i = 0; i < nN - 1; i++) {
		for (unsigned int j = i + 1; j < nN; j++) {
			dist = pbc_length_ji2(j, i);
			if (dist > cut2)
				continue;
			//cout << i << '\t' << Xp[i] << '\t' << j << '\t' << Xp[j] << '\t'					<< dist << endl;
			r = sqrt(dist);
			r_1 = 1.0 / r;
			r_2 = r_1 * r_1;
			r_6 = r_2 * r_2 * r_2;
			r_12 = r_6 * r_6;
			switchVal = 1, dSwitchVal = 0;
			if (dist > switch2) {
				c2 = cut2 - dist;
				c4 = c2 * (cut2 + 2 * dist - 3.0 * switch2);
				switchVal = c2 * c4 * c1;
				dSwitchVal = c3 * r * (c2 * c2 - c4);
				//cout << "d " << switchVal << '\t' << dSwitchVal << endl;
			}

			AmBterm = (vdwA * r_6 - vdwB) * r_6;
			Evdw += switchVal * AmBterm;
			force_r = (switchVal * 6.0 * (vdwA * r_12 + AmBterm) * r_1
					- AmBterm * dSwitchVal) * r_1;

			addF_ji(j, i, force_r);
		}
	} // Loop over self atoms

	Arg_io[AI_HAMI] = Evdw;
}

int LqFluid::Delt() {
	double force_r, AmBterm, r, r_1, r_2, r_6, r_12, switchVal, dSwitchVal, c1,
			c3, c2, c4, dist;

	//nCurTrial = ran.Number(0, nN - 1);

	//save old cor.
	CorVar[0] = Xp[nCurTrial];
	CorVar[1] = Yp[nCurTrial];
	//save old force
	FVar[0] = FXp[nCurTrial];
	FVar[1] = FYp[nCurTrial];
	//reset force
	//cout << FXp[nCurTrial] << 'v' << FYp[nCurTrial] << endl;
	FXp[nCurTrial] = 0;
	FYp[nCurTrial] = 0;

	DEvdw[0] = 0;
	c1 = 1.0 / (cut2 - switch2);
	c1 = c1 * c1 * c1;
	c3 = 4 * c1;

	for (unsigned int i = 0; i < nN; i++) {
		if (i != nCurTrial) {
			dist = pbc_length_ji2(i, nCurTrial);
			if (dist > cut2)
				continue;
			//cout << i << '\t' << j << '\t' << dist << endl;
			r = sqrt(dist);
			r_1 = 1.0 / r;
			r_2 = r_1 * r_1;
			r_6 = r_2 * r_2 * r_2;
			r_12 = r_6 * r_6;
			switchVal = 1, dSwitchVal = 0;
			if (dist > switch2) {
				c2 = cut2 - dist;
				c4 = c2 * (cut2 + 2 * dist - 3.0 * switch2);
				switchVal = c2 * c4 * c1;
				dSwitchVal = c3 * r * (c2 * c2 - c4);
				//cout << "d " << switchVal << '\t' << dSwitchVal << endl;
			}

			AmBterm = (vdwA * r_6 - vdwB) * r_6;
			DEvdw[0] += switchVal * AmBterm;
			force_r = (switchVal * 6.0 * (vdwA * r_12 + AmBterm) * r_1
					- AmBterm * dSwitchVal) * r_1;
			//cout << force_r << '\t';
			addF_ji(i, nCurTrial, force_r);
			//cout << FXp[nCurTrial] << '=' << FYp[nCurTrial] << endl;
		}
	} // Loop over self atoms
	//cout << endl;
	//cout << FXp[nCurTrial] << '^' << FYp[nCurTrial] << endl;

	Xp[nCurTrial] += ran.Real(-1, 1) * dMax;
	Yp[nCurTrial] += ran.Real(-1, 1) * dMax;
	//periodic boundary condition
	Xp[nCurTrial] = Xp[nCurTrial] - pbcx * nint(Xp[nCurTrial], pbcx);
	Yp[nCurTrial] = Yp[nCurTrial] - pbcy * nint(Yp[nCurTrial], pbcy);

	FXp[nCurTrial] = 0;
	FYp[nCurTrial] = 0;
	DEvdw[1] = 0;
	for (unsigned int i = 0; i < nN; i++) {
		if (i != nCurTrial) {
			dist = pbc_length_ji2(i, nCurTrial);
			if (dist > cut2)
				continue;
			//cout << i << '\t' << j << '\t' << dist << endl;
			r = sqrt(dist);
			r_1 = 1.0 / r;
			r_2 = r_1 * r_1;
			r_6 = r_2 * r_2 * r_2;
			r_12 = r_6 * r_6;
			switchVal = 1, dSwitchVal = 0;
			if (dist > switch2) {
				c2 = cut2 - dist;
				c4 = c2 * (cut2 + 2 * dist - 3.0 * switch2);
				switchVal = c2 * c4 * c1;
				dSwitchVal = c3 * r * (c2 * c2 - c4);
				//cout << "d " << switchVal << '\t' << dSwitchVal << endl;
			}

			AmBterm = (vdwA * r_6 - vdwB) * r_6;
			DEvdw[1] += switchVal * AmBterm;
			force_r = (switchVal * 6.0 * (vdwA * r_12 + AmBterm) * r_1
					- AmBterm * dSwitchVal) * r_1;
			//cout << force_r << '\t';
			addF_ji(i, nCurTrial, force_r);

		}
	} // Loop over self atoms
	//cout << FXp[nCurTrial] << '=' << FYp[nCurTrial] << endl;
	Arg_io[AI_DELT_HAMI] = DEvdw[1] - DEvdw[0];//new-old;
	return nCurTrial;
}

void LqFluid::MetroplisSweepBeta(int sweeps) {

	for (int CurS = 0; CurS < sweeps; CurS++) {
		for (nCurTrial = 0; nCurTrial < nN; nCurTrial++) {
			Delt();
			if (Arg_io[AI_DELT_HAMI] > 0) {
				//cout << delEvdw << '\t' << Arg_io[AI_RANDOM] << '\t'			<< Arg_io[AI_EXP_DELT_HAMI] << '\t';
				if (ran.Real() > exp(-1 * Arg_io[AI_DELT_HAMI] * dBeta)) {
					Xp[nCurTrial] = CorVar[0];//be refused and recover it
					Yp[nCurTrial] = CorVar[1];
					FXp[nCurTrial] = FVar[0];
					FYp[nCurTrial] = FVar[1];
				} //< accepted
				else {
					Arg_io[AI_HAMI] += Arg_io[AI_DELT_HAMI];
					uAccepted++;
				}
			} else {
				Arg_io[AI_HAMI] += Arg_io[AI_DELT_HAMI];
				uAccepted++;
			}
			//cout << Arg_io[AI_HAMI] << endl;
			uCurMove++;
		}
	}
	dAcceptRatio = uAccepted;
	dAcceptRatio /= uCurMove;
	if (dAcceptRatio > 0.5) {
		dMax *= 1.05;
	} else {
		dMax *= 0.95;
	}
}
void LqFluid::MetroplisSweepGCE(int sweeps) {
	double dDeltE_old2new;
	for (int CurS = 0; CurS < sweeps; CurS++) {
		for (nCurTrial = 0; nCurTrial < nN; nCurTrial++) {
			Delt();
			if (Arg_io[AI_DELT_HAMI] > 0) {
				//cout << delEvdw << '\t' << Arg_io[AI_RANDOM] << '\t'			<< Arg_io[AI_EXP_DELT_HAMI] << '\t';
				dDeltE_old2new = -Arg_io[AI_DELT_HAMI] * dBeta - dAlpha / nN
						* ((Arg_io[AI_HAMI] + Arg_io[AI_DELT_HAMI] / 2.0 - dUn
								* nN) * Arg_io[AI_DELT_HAMI]);
				if (ran.Real() > exp(dDeltE_old2new)) {
					Xp[nCurTrial] = CorVar[0];//be refused and recover it
					Yp[nCurTrial] = CorVar[1];
					FXp[nCurTrial] = FVar[0];
					FYp[nCurTrial] = FVar[1];
				} //< accepted
				else {
					Arg_io[AI_HAMI] += Arg_io[AI_DELT_HAMI];
					uAccepted++;
				}
			} else {
				Arg_io[AI_HAMI] += Arg_io[AI_DELT_HAMI];
				uAccepted++;
			}
			//cout << Arg_io[AI_HAMI] << endl;
			uCurMove++;
		}
	}
	dAcceptRatio = uAccepted;
	dAcceptRatio /= uCurMove;
	if (dAcceptRatio > 0.5) {
		dMax *= 1.05;
	} else {
		dMax *= 0.95;
	}
}

