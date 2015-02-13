/*
 * Potts.cpp
 *  Basic C++ Class for Potts model
 *  Created on: May 20, 2009
 *
 * remain precision:
 * double a=atof("-1.42343454566");
 * double b;
 * sscanf("-1.42343454566","%lf",&b);
 *
 */
#include "Potts.h"
#include "nvector.h"

#include <iostream> //< for cout
#include <cmath>
#include <ctype.h> //isdigit
using namespace std;
Potts::Potts() {

}

Potts::~Potts() {
}

void Potts::InitParameters(char q) {
	uCurMove = 0;
	cQ = q;
	States_p = NULL;
	for (int i = 0; i < AI_NUM; i++)
		Arg_io[i] = 0;
	Arg_io[AI_TEMPER] = 1;
	nGCEWF_AvgC = 10000;
	dGCEWF_SegAccumE = 0;
	nSegAccumCS = 0;
	dAccumCS = 0;
}
void Potts::SetConf(char * State_p, int dim0, int dim1) {
	nDim0 = dim0;
	nDim1 = dim1;
	nN = nDim0 * nDim1;
	nEn = nN + nN + 1;
	States_p = State_p;//delete States_p before assign new one.
}
void Potts::InitConf(int d) {
	//TODO: set random
	//ran.Srand(123487);
	if (States_p == NULL) {
		cout << "#Error: Set configures first by \'SetConf()\'" << endl;
		return;
	}
	nInitconf = d;
	if (d == 0) {
		for (unsigned int n = 0; n < nN; n++) {
			States_p[n] = ran.Number(1, cQ);
		}
	} else if (d > 0) {//TODO:NO guaranty to get anti-ground state:"four-colour problem"
		int q1, q0;
		for (q1 = cQ; q1 > 0; q1--)
			if (nDim1 % q1 != 1)
				break;
		for (q0 = cQ; q0 > 0; q0--)
			if (nDim0 % q0 != 1)
				break;

		int x0 = 0, x1 = 0;
		for (unsigned int n = 0; n < nN; n++) {
			if (n % nDim1 == 0) {
				x0++;
				x1 = 0;
			}
			States_p[n] = 1 + (x0 % q0 + x1 % q1) % cQ;
			x1++;
		}
	} else {
		for (unsigned int n = 0; n < nN; n++) {
			States_p[n] = 1;
		}
	}
}
bool Potts::LoadConf(const char * ConfName) {
	if (States_p == NULL || ConfName == NULL)
		return false;
	nvector<char> cons(States_p, nN);
	cons.load_bin(ConfName, nN, 0, false);//set last conf
	return cons.len == nN;
}
bool Potts::FromConf(const char * ConfName) {
	int clen = strlen(ConfName);
	int si = 0;
	int Q;
	while (isdigit(ConfName[si]) == 0 && si < clen) {
		si++;
	}
	if ((si < clen) && sscanf(ConfName + si, "%d_%d_%d~", &nDim0, &nDim1, &Q)
			== 3) {
		nInitconf = 2;
		InitParameters(Q);
		States_p = new char[nDim0 * nDim1 + 1];
		SetConf(States_p, nDim0, nDim1);
		return LoadConf(ConfName);
	} else {
		cout << "Unknow Conf. file:" << ConfName << endl;
		return false;
	}
}
bool Potts::DumpConf(const char * ConfName) {
	if (States_p == NULL || ConfName == NULL)
		return false;
	nvector<char> cons(States_p, nN);
	cons.dump_bin(ConfName);
	return true;
}

void Potts::SetArg(int site, double df) {
	if (site >= 0 && site < AI_NUM)
		Arg_io[site] = df;
}
double Potts::GetArg(int site) {
	if (site >= 0 && site < AI_NUM)
		return Arg_io[site];
	else {
		cout << "#Boundary error:" << site << endl;
		return 0;
	}
}

void Potts::Gauge() {
	nCurHami = 0;
	Arg_io[AI_MAG] = 0;
	unsigned int neib;
	for (unsigned int n = 0; n < nN; n++) {

		neib = n + 1;//right
		if (neib % nDim1 == 0)
			neib -= nDim1;
		if (States_p[n] == States_p[neib])
			nCurHami += 1;

		neib = n + nDim1;//down
		if (neib >= nN)
			neib -= nN;
		if (States_p[n] == States_p[neib])
			nCurHami += 1;
		Arg_io[AI_MAG] += States_p[n];

	}
	nCurHami *= -1;
	Arg_io[AI_HAMI] = nCurHami;
}

int Potts::Delt() {
	nDeltaE = 0;
	nCurTrial = ran.Number(0, nN - 1);
	char TrialValue = States_p[nCurTrial];

	//QVar = ran.Number(1, cQ);
	//choose new spin from one of the Q-1 other states
	QVar = TrialValue + ran.Number(1, cQ - 1);//[1,Q-1]
	if (QVar > cQ)
		QVar -= cQ; // differs from the present(old) spin state
	//TrialValue different from QVar
	unsigned int neib = nCurTrial + 1;//right
	if (neib % nDim1 == 0)
		neib -= nDim1;
	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	else if (States_p[neib] == QVar)
		nDeltaE -= 1;

	neib = nCurTrial - 1;//left
	if (nCurTrial % nDim1 == 0)
		neib += nDim1;
	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	else if (States_p[neib] == QVar)
		nDeltaE -= 1;

	neib = nCurTrial - nDim1;//up
	if (nCurTrial < nDim1)
		neib += nN;
	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	else if (States_p[neib] == QVar)
		nDeltaE -= 1;

	neib = nCurTrial + nDim1;//down
	if (neib >= nN)
		neib -= nN;
	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	else if (States_p[neib] == QVar)
		nDeltaE -= 1;

	Arg_io[AI_DELT_MAG] = QVar - States_p[nCurTrial];
	Arg_io[AI_DELT_HAMI] = nDeltaE;
	return nCurTrial;
}
void Potts::DeltS(unsigned int CTrial) {
	nDeltaE = 0;
	char TrialValue = States_p[CTrial];
	//TrialValue may be same as QVar
	unsigned int neib = CTrial + 1;//right

	if (neib % nDim1 == 0)
		neib -= nDim1;
	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	if (States_p[neib] == QVar)
		nDeltaE -= 1;

	neib = CTrial - 1;//left
	if (CTrial % nDim1 == 0)
		neib += nDim1;
	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	if (States_p[neib] == QVar)
		nDeltaE -= 1;

	neib = CTrial - nDim1;//up
	if (CTrial < nDim1)
		neib += nN;
	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	if (States_p[neib] == QVar)
		nDeltaE -= 1;

	neib = CTrial + nDim1;//down
	if (neib >= nN)
		neib -= nN;
	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	if (States_p[neib] == QVar)
		nDeltaE -= 1;

	//Arg_io[AI_DELT_HAMI] = nDeltaE;
	//return nDeltaE;
}
double Potts::QDistri(unsigned int * QDis_p, int M)//M=Q+1,Nmax return Nvec2;
{
	if (M == int(cQ) + 1) {
		for (int q = 0; q < cQ + 1; q++) {
			QDis_p[q] = 0;
		}
		for (unsigned int n = 0; n < nN; n++) {
			QDis_p[States_p[n]]++;
		}
		int max = 1;
		double N2 = QDis_p[1] * QDis_p[1];
		for (int q = 2; q < cQ + 1; q++) {
			N2 += QDis_p[q] * QDis_p[q];
			if (QDis_p[q] > QDis_p[max])
				max = q;
		}
		QDis_p[0] = QDis_p[max];//max,1,2...Q
		return sqrt(N2) / nN;//REF: Y.Ozeki et al./Physica A 321 (2003) 271-279
	} else {
		return -1;
	}
}
void Potts::SetTemperature(double T) {
	Arg_io[AI_TEMPER] = T;
	for (int dE = -4; dE <= 4; dE++) {
		_dDeltaExp[dE + 4] = exp(-1 * dE / Arg_io[AI_TEMPER]);//exp_de
	}
}
void Potts::MetroplisTrial(int Moves) {
	int CurStep = 0;
	while (CurStep < Moves) {
		Delt();
		//1 absolutely accepted,0 modestly accepted,-1  refused
		if (Arg_io[AI_DELT_HAMI] > 0) {
			Arg_io[AI_EXP_DELT_HAMI] = _dDeltaExp[nDeltaE + 4];
			//exp(-1 * Arg_io[AI_DELT_HAMI]/ Arg_io[AI_TEMPER]);//time-consuming instead of hash index
			Arg_io[AI_RANDOM] = ran.Real();//< random double value in the interval [0,1]
			if (Arg_io[AI_RANDOM] < Arg_io[AI_EXP_DELT_HAMI]) {
				States_p[nCurTrial] = QVar;//< accepted
				nCurHami += nDeltaE;
				//Arg_io[AI_MAG] += Arg_io[AI_DELT_MAG];
			} // else refused
		} else {
			States_p[nCurTrial] = QVar;//< accepted
			nCurHami += nDeltaE;
			//Arg_io[AI_MAG] += Arg_io[AI_DELT_MAG];
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
void Potts::MetroplisTrialBT(int Moves) {
	int CurStep = 0;
	while (CurStep < Moves) {
		Delt();
		if (Arg_io[AI_DELT_HAMI] > 0) {
			if (ran.Real() < exp(-1 * Arg_io[AI_DELT_HAMI] * Arg_io[AI_TEMPER])) {
				States_p[nCurTrial] = QVar;//< accepted
				nCurHami += nDeltaE;
			} //< refused

		} else {
			States_p[nCurTrial] = QVar;//< accepted
			nCurHami += nDeltaE;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
void Potts::MetroplisTrialGCE(int Moves) {
	int CurStep = 0;
	double dDeltE_old2new;

	while (CurStep < Moves) {
		Delt();
		if (Arg_io[AI_DELT_HAMI] > 0) {
			//dDeltE_old2new = _Se_p[nEn - 1 + nCurHami] - _Se_p[nEn - 1 + nCurHami + nDeltaE];
			dDeltE_old2new = -nDeltaE * dBeta - dAlpha / nN * ((2 * nCurHami
					+ nDeltaE - 2 * dUn * nN) * nDeltaE) / 2.0;
			if (ran.Real() < exp(dDeltE_old2new)) {
				States_p[nCurTrial] = QVar;//< accepted
				nCurHami += nDeltaE;
			} //else refused

		} else {
			States_p[nCurTrial] = QVar;//< accepted
			nCurHami += nDeltaE;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
//========== NVE ===========//
bool Potts::SetRefU(double u) {
	Gauge();
	if (u * nN > Arg_io[AI_HAMI]) {//starting with state with smaller energy
		dRefU = u * nN;
		return true;
	} else
		return false;
}
void Potts::MetroplisTrialNVE(int Moves) {
	int CurStep = 0;
	double NewHami;
	double ExpF = (nN - 2) / 2.0;
	while (CurStep < Moves) {
		Delt();
		//Accepted();
		NewHami = Arg_io[AI_HAMI] + Arg_io[AI_DELT_HAMI];
		if (dRefU >= NewHami) {
			if (Arg_io[AI_DELT_HAMI] > 0) {
				Arg_io[AI_EXP_DELT_HAMI] = pow((dRefU - NewHami) / (dRefU
						- Arg_io[AI_HAMI]), ExpF);
				if (ran.Real() < Arg_io[AI_EXP_DELT_HAMI]) {
					States_p[nCurTrial] = QVar;//< accepted
					Arg_io[AI_HAMI] += Arg_io[AI_DELT_HAMI];
				}//< refused
			} else {
				States_p[nCurTrial] = QVar;//< accepted
				Arg_io[AI_HAMI] += Arg_io[AI_DELT_HAMI];
			}
		} //else refused;
		uCurMove++;
		CurStep += 1;
	}
}

void Potts::MetroplisTrialCreutzNVE(int Moves) {
	int CurStep = 0;

	while (CurStep < Moves) {
		Delt();
		if (nDeltaE > 0) {
			if (nEDemon >= nDeltaE) {//allow dRefU=0;
				States_p[nCurTrial] = QVar;//< accepted
				nEDemon -= nDeltaE;
				nCurHami += nDeltaE;
			} //else refused
		} else {
			States_p[nCurTrial] = QVar;//< accepted
			nEDemon -= nDeltaE;
			nCurHami += nDeltaE;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
//========== HEAT Bath ===========//
bool Potts::SetHeatPath(double T, int *DEi, int n) {
	if (n != cQ)
		return false;
	Arg_io[AI_TEMPER] = T;
	double sde = 0;
	for (int mE = 0; mE <= 4; mE++) {//Q>3
		_dDeltaExp[mE] = exp(-1 * mE / Arg_io[AI_TEMPER]);//exp_de
		sde += _dDeltaExp[mE];
	}
	for (int mE = 0; mE <= 4; mE++) {//normalization
		_dDeltaExp[mE] /= sde;
		cout << _dDeltaExp[mE] << '\t';
	}
	for (int mE = 1; mE <= 4; mE++) {//Accum
		_dDeltaExp[mE] += _dDeltaExp[mE - 1];
	}
	QDeltaEi = DEi;
	return true;
}
int Potts::HeatPathDelt() {
	nCurTrial = ran.Number(0, nN - 1);
	//char TrialValue = States_p[nCurTrial];
	unsigned int Dir[4];
	Dir[0] = nCurTrial + 1;//right
	if (Dir[0] % nDim1 == 0)
		Dir[0] -= nDim1;

	Dir[1] = nCurTrial - 1;//left
	if (nCurTrial % nDim1 == 0)
		Dir[1] += nDim1;

	Dir[2] = nCurTrial - nDim1;//up
	if (nCurTrial < nDim1)
		Dir[2] += nN;

	Dir[3] = nCurTrial + nDim1;//down
	if (Dir[3] >= nN)
		Dir[3] -= nN;

	for (QVar = 1; QVar <= cQ; QVar++) {
		nDeltaE = 0;
		for (int i = 0; i < 4; i++)
			if (States_p[Dir[i]] == QVar)
				nDeltaE += 1;
		QDeltaEi[QVar - 1] = nDeltaE;
		cout << QDeltaEi[QVar - 1] << '\t';
	}
	cout << endl;
	return nCurTrial;
}
void Potts::HeatPathTrial(int Moves)//NVT
{

	/*
	 * int CurStep = 0;
	 while (CurStep < Moves) {
	 HeatPathDelt();
	 Arg_io[AI_RANDOM] = ran.Real();//< random double value in the interval [0,1]
	 double AccPro[cQ], sde = 0;

	 for (QVar = 1; QVar <= cQ; QVar++) {//Q>3
	 AccPro[QVar - 1] = _dDeltaExp[QDeltaEi[QVar - 1]];
	 sde += AccPro[QVar - 1];
	 }
	 for (QVar = 1; QVar <= cQ; QVar++) {//normalization,accum
	 AccPro[QVar - 1] /= sde;
	 }
	 cout << AccPro[0] << '\t';
	 for (QVar = 2; QVar <= cQ; QVar++) {
	 AccPro[QVar - 1] += AccPro[QVar - 2];
	 cout << AccPro[QVar - 1] << '\t';
	 }
	 cout << endl;

	 for (QVar = 1; QVar <= cQ; QVar++) {
	 if (Arg_io[AI_RANDOM] < AccPro[QVar - 1]) {
	 break;
	 }

	 }
	 cout << Arg_io[AI_RANDOM] << '\t' << int(QVar) << '\t' << nDeltaE
	 << endl;
	 Arg_io[AI_EXP_DELT_HAMI] = (nDeltaE > 0) ? _dDeltaExp[nDeltaE]
	 - _dDeltaExp[nDeltaE - 1] : _dDeltaExp[nDeltaE];
	 Arg_io[AI_DELT_HAMI] = QDeltaEi[States_p[nCurTrial] - 1]
	 - QDeltaEi[QVar - 1];

	 Arg_io[AI_HAMI] += Arg_io[AI_DELT_HAMI];
	 //Arg_io[AI_MAG] += Arg_io[AI_DELT_MAG];
	 States_p[nCurTrial] = QVar;//< accepted
	 uCurMove++;
	 CurStep += 1;
	 }*/
}

//========== FAFE ===========//
bool Potts::SetSefun(double * Se_p, int size) {
	if (nEn > size) {
		cerr << "#Error: Size too small:" << size << endl;
		return false;
	} else {
		_Se_p = Se_p;
		return true;
	}
}
void Potts::MetroplisTrialSe(int Moves) {
	int CurStep = 0;
	double dDeltE_old2new;
	while (CurStep < Moves) {
		Delt();
		if (Arg_io[AI_DELT_HAMI] > 0) {
			dDeltE_old2new = _Se_p[nEn - 1 + nCurHami] - _Se_p[nEn - 1
					+ nCurHami + nDeltaE];
			if (ran.Real() < exp(dDeltE_old2new)) {
				States_p[nCurTrial] = QVar;//< accepted
				nCurHami += nDeltaE;
			} //else refused
		} else {
			States_p[nCurTrial] = QVar;//< accepted
			nCurHami += nDeltaE;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
//////////////////////////////////
bool Potts::SetWefun(double * We_p, int size) {
	if (nEn > size) {
		cerr << "#Error: Size too small:" << size << endl;
		return false;
	} else {
		_We_p = We_p;
		return true;
	}
}
void Potts::MetroplisTrialWe(int Moves) {
	int CurStep = 0;

	while (CurStep < Moves) {
		Delt();
		if (Arg_io[AI_DELT_HAMI] > 0) {
			Arg_io[AI_EXP_DELT_HAMI] = _We_p[nEn - 1 + nCurHami + nDeltaE]
					/ _We_p[nEn - 1 + nCurHami];
			Arg_io[AI_RANDOM] = ran.Real();//< random double value in the interval [0,1]
			if (Arg_io[AI_RANDOM] < Arg_io[AI_EXP_DELT_HAMI]) {
				States_p[nCurTrial] = QVar;//< accepted
				nCurHami += nDeltaE;
			} //else refused
		} else {
			States_p[nCurTrial] = QVar;//< accepted
			nCurHami += nDeltaE;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
//======== Wang-Landau ==========//

void Potts::SetWL(double * lnGe_p, unsigned long *geN_p, int size,
		double lnfPrecision, double flatRate) {
	if (nEn > size) {
		cout << "#Error: Size too small:" << size << endl;
		return;
	}
	LnfPrecision = lnfPrecision;
	FlatRate = flatRate;
	LnGe_p = lnGe_p;
	GeN_p = geN_p;
	Lnf = 1;
	for (int n = 0; n < nEn; n++) {
		GeN_p[n] = 0;
		LnGe_p[n] = 0;
	}
}

bool Potts::WLCheckBelowAVG() {
	double dAvg = uCurMove;
	/*
	 for (int n = 6; n < nEn; n++) {
	 dAvg += GeN_p[n];
	 }
	 dAvg += GeN_p[4];
	 dAvg += GeN_p[0];
	 */
	dAvg /= nN + nN - 3;
	dAvg *= FlatRate;
	if (GeN_p[0] < dAvg)
		return false;
	if (GeN_p[4] < dAvg)
		return false;
	for (int n = 6; n < nEn; n++) {
		if (GeN_p[n] < dAvg)
			return false;
	}

	return true;
}
bool Potts::WLCheckBelowAVGBound() {
	double dAvg = uCurMove;
	/*
	 for (int n = 6; n < nEn; n++) {
	 dAvg += GeN_p[n];
	 }
	 dAvg += GeN_p[4];
	 dAvg += GeN_p[0];
	 */
	/*
	 * -2N,_,_,_,-2N-4,_,-2N-5,2N-6,..,0
	 * 0,...4,,6
	 */
	int lEi = nEn - 1 + nLowEBound;
	int hEi = nEn - 1 + nHighEBound;
	//assume (lEi >= 6)
	dAvg /= nHighEBound - nLowEBound + 1;
	dAvg *= FlatRate;
	for (int n = lEi; n <= hEi; n++) {
		if (GeN_p[n] < dAvg)
			return false;
	}
	return true;
}
void Potts::WLResetGeN() {
	uCurMove = 0;

	for (int n = 1; n < nEn; n++) {
		GeN_p[n] = 0;
		LnGe_p[n] -= LnGe_p[0];//save precisions as far as possible
	}
	GeN_p[0] = 0;
	LnGe_p[5] = 0;
	for (int n = 0; n < 4; n++)
		LnGe_p[n] = 0;

}

void Potts::WangLandauTrial(int Moves) {
	int CurStep = 0;
	int E1, E2;
	double dDeltLnGE_old2new;

	while (CurStep < Moves) {
		Delt();
		E1 = nEn - 1 + nCurHami;
		E2 = E1 + nDeltaE;
		dDeltLnGE_old2new = LnGe_p[E1] - LnGe_p[E2];
		if (dDeltLnGE_old2new < 0) {
			Arg_io[AI_RANDOM] = ran.Real();//
			Arg_io[AI_EXP_DELT_HAMI] = exp(dDeltLnGE_old2new);
			if (Arg_io[AI_RANDOM] < Arg_io[AI_EXP_DELT_HAMI]) {
				LnGe_p[E2] += Lnf;
				GeN_p[E2] += 1;
				States_p[nCurTrial] = QVar;//< accepted
				nCurHami += nDeltaE;
				//nAccepted = 0;
			} else {
				LnGe_p[E1] += Lnf;
				GeN_p[E1] += 1;
				//nAccepted = -1;
			}
		} else {
			LnGe_p[E2] += Lnf;
			GeN_p[E2] += 1;
			States_p[nCurTrial] = QVar;//< accepted
			nCurHami += nDeltaE;
			//nAccepted = 1;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
void Potts::WangLandauTrialBound(int Moves) {
	int CurStep = 0;
	int E1, E2;
	double dDeltLnGE_old2new;

	while (CurStep < Moves) {
		Delt();//Assume nLowEBound<=nCurHami<=nHighEBound
		if ((nCurHami + nDeltaE) < nLowEBound || (nCurHami + nDeltaE)
				> nHighEBound)
			continue;

		E1 = nEn - 1 + nCurHami;
		E2 = E1 + nDeltaE;
		dDeltLnGE_old2new = LnGe_p[E1] - LnGe_p[E2];

		if (dDeltLnGE_old2new < 0) {
			Arg_io[AI_RANDOM] = ran.Real();//
			Arg_io[AI_EXP_DELT_HAMI] = exp(dDeltLnGE_old2new);
			if (Arg_io[AI_RANDOM] < Arg_io[AI_EXP_DELT_HAMI]) {
				LnGe_p[E2] += Lnf;
				GeN_p[E2] += 1;
				States_p[nCurTrial] = QVar;//< accepted
				nCurHami += nDeltaE;
				//nAccepted = 0;
			} else {
				LnGe_p[E1] += Lnf;
				GeN_p[E1] += 1;
				//nAccepted = -1;
			}
		} else {
			LnGe_p[E2] += Lnf;
			GeN_p[E2] += 1;
			States_p[nCurTrial] = QVar;//< accepted
			nCurHami += nDeltaE;
			//nAccepted = 1;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
void Potts::WarmupToBound() {
	InitConf(-1);//warm up
	Gauge();
	while (true) {
		Delt();
		States_p[nCurTrial] = QVar;//< accepted,T=0
		nCurHami += nDeltaE;
		if ((nLowEBound <= nCurHami) && (nCurHami <= nHighEBound)) {
			break;
		}
	}
}
//======== Multicanonical ==========//

int Potts::SetMul(double * beta_p, double *se_p, double *gn_p,
		unsigned long *hist_p, int size, double h)//len of beta_p,_Se_p,Hist_p,Gn_p =nEn,
{
	if (beta_p == NULL || (gn_p == NULL) || (hist_p == NULL) || (se_p == NULL))
		return -1;

	if (size != nEn)
		return -2;
	if (h > 1 || h < 0)
		return -3;
	Beta_p = beta_p;
	_Se_p = se_p;
	Gn_p = gn_p;
	Hist_p = hist_p;
	dh = h;
	return 1;
}
void Potts::MulReSet() {
	int bLei = nEn - 1 + nLowEBound;
	int bHei = nEn - 1 + nHighEBound;

	_Se_p[bLei] = 0;
	Hist_p[bLei] = 0;
	for (int ei = bLei + 1; ei <= bHei; ei++) {
		_Se_p[ei] = _Se_p[ei - 1] + Beta_p[ei - 1];
		Hist_p[ei] = 0;
	}

}
void Potts::MetroplisTrialSeBound(int Moves)//[nLowEBound, nHighEBound]
{
	int CurStep = 0;
	double dDeltE_old2new;

	while (CurStep < Moves) {
		Delt();//Assume nLowEBound<=nCurHami<=nHighEBound
		if ((nCurHami + nDeltaE) < nLowEBound || (nCurHami + nDeltaE)
				> nHighEBound)
			continue;
		if (nDeltaE > 0) {
			dDeltE_old2new = _Se_p[nEn - 1 + nCurHami] - _Se_p[nEn - 1
					+ nCurHami + nDeltaE];
			if (ran.Real() < exp(dDeltE_old2new)) {
				States_p[nCurTrial] = QVar;//< accepted
				nCurHami += nDeltaE;
			} //else refused
		} else {
			States_p[nCurTrial] = QVar;//< accepted
			nCurHami += nDeltaE;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
bool Potts::MulReweight() {
	double gn0, cap_gn0;
	int bLei = nEn - 1 + nLowEBound;
	int bHei = nEn - 1 + nHighEBound;

	for (int ei = bLei; ei < bHei; ei++) {
		//cout << Hist_p[ei + 1] << '\t' << Hist_p[ei] << endl;
		gn0 = (Hist_p[ei + 1] + Hist_p[ei]);

		if (gn0 != 0)
			gn0 = Hist_p[ei + 1] * Hist_p[ei] / gn0;

		Gn_p[ei] += gn0;
		cap_gn0 = gn0 / Gn_p[ei];
		if (cap_gn0 > 0) {
			Beta_p[ei] += cap_gn0 * log(((Hist_p[ei + 1] > 0) ? Hist_p[ei + 1]
					: dh) / ((Hist_p[ei] > 0) ? Hist_p[ei] : dh));
		}
	}
	bool allNonZero = true;
	for (int ei = bLei; ei < bHei; ei++) {
		if (Hist_p[ei] <= nFlatNonzero) {
			allNonZero = false;
			break;
		}
	}
	return allNonZero;

}
//=========================//

void Potts::WolffCluster(int Moves) {
	int CurStep = 0;
	int sp, oldspin, newspin;
	unsigned int neib;
	//int oldbonds = 0, newbonds = 0;
	while (CurStep < Moves) {
		nCurTrial = ran.Number(0, nN - 1);
		ClusterStack_p[0] = nCurTrial;
		sp = 1;
		oldspin = States_p[nCurTrial];
		// chose a new spin state
		newspin = oldspin + ran.Number(1, cQ - 1);//[1,Q-1]
		if (newspin > cQ)
			newspin -= cQ; // differs from the present(old) spin state

		States_p[nCurTrial] = newspin;
		while (sp) {
			nCurTrial = ClusterStack_p[--sp];
			neib = nCurTrial + 1;//right
			if (neib % nDim1 == 0)
				neib -= nDim1;
			if (States_p[neib] == oldspin) {
				if (ran.Real() < dPadd) {
					ClusterStack_p[sp++] = neib;
					States_p[neib] = newspin;
				}//oldbonds++;
			}
			neib = nCurTrial - 1;//left
			if (nCurTrial % nDim1 == 0)
				neib += nDim1;
			if (States_p[neib] == oldspin) {
				if (ran.Real() < dPadd) {
					ClusterStack_p[sp++] = neib;
					States_p[neib] = newspin;
				}
			}
			neib = nCurTrial - nDim1;//up
			if (nCurTrial < nDim1)
				neib += nN;
			if (States_p[neib] == oldspin) {
				if (ran.Real() < dPadd) {
					ClusterStack_p[sp++] = neib;
					States_p[neib] = newspin;
				}
			}

			neib = nCurTrial + nDim1;//down
			if (neib >= nN)
				neib -= nN;
			if (States_p[neib] == oldspin) {
				if (ran.Real() < dPadd) {
					ClusterStack_p[sp++] = neib;
					States_p[neib] = newspin;
				}

			}
		}
		uCurMove++;
		CurStep += 1;
	}
}
void Potts::ClusterDelt(unsigned int CTrial) {//QVar=newspin
//nDeltaE = 0;
	char TrialValue = States_p[CTrial];
	//TrialValue different from QVar
	unsigned int neib = CTrial + 1;//right
	if (neib % nDim1 == 0)
		neib -= nDim1;

	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	else if (States_p[neib] == QVar)
		nDeltaE -= 1;

	neib = CTrial - 1;//left
	if (CTrial % nDim1 == 0)
		neib += nDim1;

	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	else if (States_p[neib] == QVar)
		nDeltaE -= 1;

	neib = CTrial - nDim1;//up
	if (CTrial < nDim1)
		neib += nN;

	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	else if (States_p[neib] == QVar)
		nDeltaE -= 1;

	neib = CTrial + nDim1;//down
	if (neib >= nN)
		neib -= nN;

	if (States_p[neib] == TrialValue)
		nDeltaE += 1;
	else if (States_p[neib] == QVar)
		nDeltaE -= 1;

	//Arg_io[AI_DELT_HAMI] = nDeltaE;
	//return nDeltaE;
}
void Potts::WolffClusterGCE(int Moves) {
	//set first dBeta, dAlpha, dUn;
	int CurStep = 0;
	int sp = -1, oldspin, newspin;

	unsigned int neib;
	//int oldbonds = 0, newbonds = 0;
	while (CurStep < Moves) {
		sp = 0;
		csp = 0;
		nDeltaE = 0;
		nCurTrial = ran.Number(0, nN - 1);
		ClusterStack_p[sp++] = nCurTrial;
		dPadd = 1 - exp(-dBeta - dAlpha * (dGCEWF_AvgE - dUn * nN) / nN);
		oldspin = States_p[nCurTrial];
		// choose a new spin state
		newspin = oldspin + ran.Number(1, cQ - 1);//[1,Q-1]
		if (newspin > cQ)
			newspin -= cQ; // differs from the present(old) spin state
		QVar = newspin;
		ClusterDelt(nCurTrial);
		ClusterSeq_p[csp++] = nCurTrial;
		States_p[nCurTrial] = newspin;
		while (sp) {
			nCurTrial = ClusterStack_p[--sp];
			neib = nCurTrial + 1;//right
			if (neib % nDim1 == 0)
				neib -= nDim1;

			if (States_p[neib] == oldspin) {
				if (ran.Real() < dPadd) {
					ClusterStack_p[sp++] = neib;
					ClusterDelt(neib);
					ClusterSeq_p[csp++] = neib;
					States_p[neib] = newspin;
				}
			}
			neib = nCurTrial - 1;//left
			if (nCurTrial % nDim1 == 0)
				neib += nDim1;
			if (States_p[neib] == oldspin) {
				if (ran.Real() < dPadd) {
					ClusterStack_p[sp++] = neib;
					ClusterDelt(neib);
					ClusterSeq_p[csp++] = neib;
					States_p[neib] = newspin;
				}
			}
			neib = nCurTrial - nDim1;//up
			if (nCurTrial < nDim1)
				neib += nN;
			if (States_p[neib] == oldspin) {
				if (ran.Real() < dPadd) {
					ClusterStack_p[sp++] = neib;
					ClusterDelt(neib);
					ClusterSeq_p[csp++] = neib;
					States_p[neib] = newspin;
				}
			}

			neib = nCurTrial + nDim1;//down
			if (neib >= nN)
				neib -= nN;
			if (States_p[neib] == oldspin) {
				if (ran.Real() < dPadd) {
					ClusterStack_p[sp++] = neib;
					ClusterDelt(neib);
					ClusterSeq_p[csp++] = neib;
					States_p[neib] = newspin;
				}

			}

		}//cluster formed now
		if (ran2.Real() < exp(-dAlpha * nDeltaE * (nCurHami + nDeltaE / 2.0
				- dGCEWF_AvgE) / nN)) {//accepted
			nCurHami += nDeltaE;
			Arg_io[AI_HAMI] = nCurHami;

		} else {//refused
			for (int p = 0; p < csp; p++) {
				States_p[ClusterSeq_p[p]] = oldspin;
			}//roll back.
		}
		//nDeltaE = 0;
		uCurMove++;
		CurStep += 1;
		dGCEWF_SegAccumE += nCurHami;
		nSegAccumCS += csp;
		if (uCurMove % nGCEWF_AvgC == 0) {
#ifdef _DEBUG
			cout << uCurMove << '\t' << csp << '\t' << dPadd << '\t' << exp(
					-dAlpha * nDeltaE
					* (nCurHami + nDeltaE / 2.0 - dGCEWF_AvgE) / nN)
			<< '\t' << nDeltaE << '\t' << dGCEWF_AvgE <<'\t'<<nSegAccumCS / double(nGCEWF_AvgC)<< endl;//2.0!!
#endif

			dGCEWF_AvgE = dGCEWF_SegAccumE / nGCEWF_AvgC;
			dGCEWF_SegAccumE = 0;

			dAccumCS += nSegAccumCS / double(nGCEWF_AvgC);
			nSegAccumCS = 0;

		}
	}
}
