/*
 * SeFun.cpp
 *
 *  Created on: Mar 9, 2009
 */

#include "SeFun.h"
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
using namespace std;
////////////////////////
#ifndef M_PI
#define M_PI 3.1415926535897931
#endif
SeFun::SeFun() {
	bSmooth = true;
}
SeFun::~SeFun() {
}
void SeFun::ReSetFun(int nNum0, int nNum1, double *dTE) {
	_TELine_l.clear();
	for (int TempIndex = 0; TempIndex < nNum0; TempIndex++) {
		TEP tep;
		tep.T = dTE[nNum1 * TempIndex];
		tep.E = dTE[nNum1 * TempIndex + 1];
		tep.Sigma = dTE[nNum1 * TempIndex + 2];
		tep.A = dTE[nNum1 * TempIndex + 3];
		_TELine_l.push_back(tep);
		//TODO: remove print
		//cout << _TELine_l[TempIndex].T << '\t' << _TELine_l[TempIndex].E<< '\t' << _TELine_l[TempIndex].Sigma << endl;
	}
}
void SeFun::InitFun(double lflatEx, double hflatEx, bool smooth) {
	vector<TEP>::iterator next, prev;
	prev = _TELine_l.begin();
	if (prev == _TELine_l.end()) {
		cout << "vector<TEP> is empty" << endl;
		return;
	}

	bSmooth = smooth;
	double beta = 0;
	dLflatEx = prev->Sigma * lflatEx;
	prev->SE = dLflatEx / prev->T;//set low-temperature end rather than high-temperature end
	next = prev;
	next++;
	while (next != _TELine_l.end()) {
		next->SE = prev->SE + _DIntegralBeta(prev, next, next->E, &beta);
		prev = next;
		next++;
	}
	next--;
	dHflatEx = next->Sigma * hflatEx;
}

double SeFun::SEFunction(double energy, double *beta, double High_t,
		double Low_t) {
	double SE_t = 0;

	vector<TEP>::iterator next, last, prev;
	next = _TELine_l.begin();
	last = _TELine_l.end();
	last--;
	//cout << "1" << next->E << '\t' << next->FlatEx << '\t' << energy << endl;
	if (energy < (next->E - dLflatEx)) {
		SE_t = (energy - next->E + dLflatEx) / High_t;
		*beta = 1.0 / High_t;

	} else if (energy <= next->E) {
		SE_t = (energy - next->E + dLflatEx) / next->T;
		*beta = 1.0 / next->T;
	} else {
		while (next != last && energy > next->E) {
			next++;
		}
		if (energy > last->E + dHflatEx) {
			SE_t = last->SE + (dHflatEx - last->E) / last->T + (energy
					- dHflatEx) / Low_t;// flat high;
			*beta = 1.0 / Low_t;
		} else if (energy > last->E) {
			SE_t = last->SE + (energy - last->E) / last->T;// flat high;
			*beta = 1.0 / last->T;
		} else {
			prev = next;
			prev--;
			SE_t = prev->SE + _DIntegralBeta(prev, next, energy, beta);

		}
	}

	return SE_t;
}
double SeFun::_DIntegralBeta(vector<TEP>::iterator begin,
		vector<TEP>::iterator end, double energy, double *beta) {
	double SE_t = 0, dSlope_t = 0;
	if (DOUBLE_EQ(begin->T, end->T)) {
		*beta = 1.0 / begin->T;
		SE_t = (energy - begin->E) / begin->T;
	} else {
		dSlope_t = (end->T - begin->T) / (end->E - begin->E);
		*beta = 1.0 / (begin->T + dSlope_t * (energy - begin->E));
		SE_t = log(abs(1.0 / *beta / begin->T)) / dSlope_t;
	}

	if (bSmooth) {
		double slope1_t = begin->A - 1.0 / (begin->Sigma * begin->Sigma);
		double slope2_t = end->A - 1.0 / (end->Sigma * end->Sigma);

		double omega = M_PI / (end->E - begin->E);
		double b1 = slope1_t + dSlope_t / begin->T / begin->T;
		double b2 = slope2_t + dSlope_t / end->T / end->T;
		double coeff1 = b1 - b2;
		coeff1 /= 2.0 * omega;
		double coeff2 = b1 + b2;
		coeff2 /= 4.0 * omega;
		double theta1 = (energy - begin->E) * omega;
		*beta += coeff1 * sin(theta1) + coeff2 * sin(2.0 * theta1);

		SE_t -= (coeff1 / omega * (cos(theta1) - 1.0)) + coeff2 / (2.0 * omega)
				* (cos(2.0 * theta1) - 1.0);
	}
	return SE_t;
}
