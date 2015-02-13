/*
 * Ising.cpp
 *
 *  Created on: May 20, 2009
 *   0---N-->
 * (0,0)----y,dim1--->
 *   |
 *   |
 *   x,
 *  dim0
 *   |
 *   v
 * Energy Levels(N+1):
 * -2N,,-2N-8,2N-12,..,2N-8,2N
 * for(int n=0;n<=N;n++)
 * 		E[n]= 4*n-2*N;
 * hisE[(2N + hami) / 4)]++
 * M:lev=N+1
 * -N,-N-2,...,0,2,...,N-2,N
 */
#include "Ising.h"
#include "nvector.h"
#include <iostream> //< for cout
#include <cmath>

using namespace std;

Ising::Ising() {
}

Ising::~Ising() {
}
void Ising::InitParameters() {
	uCurMove = 0;
	States_p = NULL;
	for (int i = 0; i < AI_NUM; i++)
		Arg_io[i] = 0;
	Arg_io[AI_TEMPER] = 1;
}
void Ising::SetConf(char * State_p, int dim0, int dim1) {
	nDim0 = dim0;
	nDim1 = dim1;
	nN = nDim0 * nDim1;
	nEn = nN + 1;
	/*  nEn/4 ->min non-zero Delta E=4
	 *2*nN+ 1 +2*nN                 nDim0,nDim0:both even
	 *2*nN+ 1 +2*nN-2*both_odd_dim  nDim0,nDim0:both odd
	 *2*nN+ 1 +2*nN-2*even_dim      nDim0,nDim0:one of them odd
	 */
	States_p = State_p;
}
void Ising::InitConf(int d) {
	//ran.Srand(123487);
	if (States_p == NULL) {
		cout << "#Error: Set configures first by \'SetConf()\'" << endl;
		return;
	}
	nInitconf = d;
	if (d == 0) {
		for (unsigned int n = 0; n < nN; n++) {
			States_p[n] = (ran.Real() > 0.5) ? 1 : -1;
		}
	} else if (d > 0) {//TODO:NO guaranty to get anti-ground state:"four-colour problem"
		int m, k = 1;
		for (unsigned int x = 0; x < nDim0; x++) {
			k *= -1;
			m = k;
			for (unsigned int y = 0; y < nDim1; y++) {
				States_p[x * nDim1 + y] = m;
				m *= -1;
			}
		}
	} else {
		for (unsigned int n = 0; n < nN; n++) {
			States_p[n] = -1;
		}
	}
}
bool Ising::FromConf(const char * ConfName) {
	int clen = strlen(ConfName);
	int si = 0;
	while (isdigit(ConfName[si]) == 0 && si < clen) {
		si++;
	}
	if ((si < clen) && sscanf(ConfName + si, "%d_%d~", &nDim0, &nDim1) == 2) {
		nInitconf = 2;
		InitParameters();
		States_p = new char[nDim0 * nDim1 + 1];
		SetConf(States_p, nDim0, nDim1);
		return LoadConf(ConfName);
	} else {
		cout << "Unknow Conf. file:" << ConfName << endl;
		return false;
	}
}
bool Ising::DumpConf(const char * ConfName) {
	if (States_p == NULL || ConfName == NULL)
		return false;
	nvector<char> cons(States_p, nN);
	cons.dump_bin(ConfName);
	return true;
}
bool Ising::LoadConf(const char * ConfName) {
	if (States_p == NULL || ConfName == NULL)
		return false;
	nvector<char> cons(States_p, nN);
	cons.load_bin(ConfName, nN, 0, false);//set last conf
	return cons.len == nN;
}

void Ising::SetArg(int site, double df) {
	if (site >= 0 && site < AI_NUM)
		Arg_io[site] = df;
}
double Ising::GetArg(int site) {
	if (site >= 0 && site < AI_NUM)
		return Arg_io[site];
	else {
		cout << "#Boundary error:" << site << endl;
		return 0;
	}
}

void Ising::Gauge() {
	nCurHami = 0;
	nCurMag = 0;
	unsigned int neib;
	for (unsigned int n = 0; n < nN; n++) {
		neib = n + 1;//right
		if (neib % nDim1 == 0)
			neib -= nDim1;
		nCurHami += States_p[n] * States_p[neib];

		neib = n + nDim1;//down
		if (neib >= nN)
			neib -= nN;
		nCurHami += States_p[n] * States_p[neib];
		nCurMag += States_p[n];
	}
	nCurHami *= -1;
	Arg_io[AI_HAMI] = nCurHami;
	Arg_io[AI_MAG] = nCurMag;
}

int Ising::Delt() {
	nDeltaE = 0;
	nCurTrial = ran.Number(0, nN - 1);

	unsigned int neib = nCurTrial + 1;//right
	if (neib % nDim1 == 0)
		neib -= nDim1;
	nDeltaE += States_p[neib];

	neib = nCurTrial - 1;//left
	if (nCurTrial % nDim1 == 0)
		neib += nDim1;
	nDeltaE += States_p[neib];

	neib = nCurTrial - nDim1;//up
	if (nCurTrial < nDim1)
		neib += nN;
	nDeltaE += States_p[neib];

	neib = nCurTrial + nDim1;//down
	if (neib >= nN)
		neib -= nN;
	nDeltaE += States_p[neib];

	nDeltaE *= (States_p[nCurTrial] << 1);
	nDeltaMag = -(States_p[nCurTrial] << 1);
	Arg_io[AI_DELT_MAG] = nDeltaMag;
	Arg_io[AI_DELT_HAMI] = nDeltaE;
	return nCurTrial;
}

void Ising::DeltS(unsigned int CTrial) {
	nDeltaE = 0;

	unsigned int neib = CTrial + 1;//right
	if (neib % nDim1 == 0)
		neib -= nDim1;
	nDeltaE += States_p[neib];

	neib = CTrial - 1;//left
	if (CTrial % nDim1 == 0)
		neib += nDim1;
	nDeltaE += States_p[neib];

	neib = CTrial - nDim1;//up
	if (CTrial < nDim1)
		neib += nN;
	nDeltaE += States_p[neib];

	neib = CTrial + nDim1;//down
	if (neib >= nN)
		neib -= nN;
	nDeltaE += States_p[neib];

	nDeltaE *= (2 * States_p[CTrial]);

	//	Arg_io[AI_DELT_HAMI] = nDeltaE;
}
void Ising::SetTemperature(double T) {
	Arg_io[AI_TEMPER] = T;//dE=2JS@s4,[-8,-4,0,4,8]
	for (int dE = -8; dE <= 8; dE += 4) {
		_dDeltaExp[2 + dE / 4] = exp(-1 * dE / Arg_io[AI_TEMPER]);//exp_de
	}
}
void Ising::MetroplisTrial(int Moves) {
	int CurStep = 0;
	while (CurStep < Moves) {
		Delt();
		//Accepted();//1 absolutely accepted,0 modestly accepted,-1  refused
		if (nDeltaE > 0) {
			Arg_io[AI_EXP_DELT_HAMI] = _dDeltaExp[2 + nDeltaE / 4];
			//exp(-1 * Arg_io[AI_DELT_HAMI]/ Arg_io[AI_TEMPER]);//time-consuming instead of hash index
			//Arg_io[AI_RANDOM] = ran.Real();//< random double value in the interval [0,1]
			if (ran.Real() < Arg_io[AI_EXP_DELT_HAMI]) {
				States_p[nCurTrial] *= -1;//< accepted
				nCurHami += nDeltaE;
				nCurMag += nDeltaMag;
				nAccepted = 0;
			} else {//< refused
				nAccepted = -1;
			}
		} else {
			States_p[nCurTrial] *= -1;//< accepted
			nCurHami += nDeltaE;
			nCurMag += nDeltaMag;
			nAccepted = 1;
		}
		uCurMove++;
		CurStep++;
	}
	Arg_io[AI_HAMI] = nCurHami;
	Arg_io[AI_MAG] = nCurMag;
}
void Ising::MetroplisTrialGCE(int Moves) {
	int CurStep = 0;
	double dDeltE_old2new;

	while (CurStep < Moves) {
		Delt();
		if (nDeltaE > 0) {
			//dDeltE_old2new = _Se_p[nEn - 1 + nCurHami] - _Se_p[nEn - 1 + nCurHami + nDeltaE];
			dDeltE_old2new = -nDeltaE * dBeta - dAlpha / nN * ((2 * nCurHami
					+ nDeltaE - 2 * dUn * nN) * nDeltaE) / 2.0;
			if (ran.Real() < exp(dDeltE_old2new)) {
				States_p[nCurTrial] *= -1;//< accepted
				nCurHami += nDeltaE;
				nCurMag += nDeltaMag;

			} //else refused

		} else {
			States_p[nCurTrial] *= -1;//< accepted
			nCurHami += nDeltaE;
			nCurMag += nDeltaMag;
		}
		uCurMove++;
		CurStep++;
	}
	Arg_io[AI_HAMI] = nCurHami;
	Arg_io[AI_MAG] = nCurMag;
}
void Ising::WolffCluster(int Moves) {
	int CurStep = 0;
	int sp, oldspin, newspin;
	unsigned int neib;
	//int oldbonds = 0, newbonds = 0;
	while (CurStep < Moves) {
		nCurTrial = ran.Number(0, nN - 1);
		ClusterStack_p[0] = nCurTrial;
		sp = 1;
		oldspin = States_p[nCurTrial];
		newspin = -oldspin;
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
		CurStep++;
	}
}
void Ising::ClusterDelt(unsigned int CTrial) {
	//nDeltaE = 0;
	unsigned int neib = CTrial + 1;//right
	int nDE = 0;
	if (neib % nDim1 == 0)
		neib -= nDim1;
	nDE += States_p[neib];

	neib = CTrial - 1;//left
	if (CTrial % nDim1 == 0)
		neib += nDim1;
	nDE += States_p[neib];

	neib = CTrial - nDim1;//up
	if (CTrial < nDim1)
		neib += nN;
	nDE += States_p[neib];

	neib = CTrial + nDim1;//down
	if (neib >= nN)
		neib -= nN;
	nDE += States_p[neib];

	nDeltaE += nDE * (States_p[CTrial] << 1);
	nDeltaMag += -(States_p[CTrial] << 1);

	//Arg_io[AI_DELT_HAMI] = nDeltaE;
	//return nDeltaE;
}
void Ising::WolffClusterGCE(int Moves) {
	//set first dBeta, dAlpha, dUn;
	int CurStep = 0;
	int sp = -1, oldspin, newspin;

	unsigned int neib;
	//int oldbonds = 0, newbonds = 0;
	while (CurStep < Moves) {
		sp = 0;
		csp = 0;
		nDeltaE = 0;
		nDeltaMag = 0;
		nCurTrial = ran.Number(0, nN - 1);
		ClusterStack_p[sp++] = nCurTrial;
		//E_v-E_u=2J(m-n),for ising
		//Ev-Eu=J(m-n),for potts

		dPadd = 1
				- exp(-2.0 * (dBeta + dAlpha * (dGCEWF_AvgE - dUn * nN) / nN));
		oldspin = States_p[nCurTrial];
		// switch to the spin state
		newspin = -oldspin;
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
			nCurMag += nDeltaMag;
			Arg_io[AI_HAMI] = nCurHami;
			Arg_io[AI_MAG] += nDeltaMag;

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
			<< '\t' << nDeltaE << '\t' << dGCEWF_AvgE << '\t'
			<< nSegAccumCS / double(nGCEWF_AvgC) << endl;//2.0!!
#endif

			dGCEWF_AvgE = dGCEWF_SegAccumE / nGCEWF_AvgC;
			dGCEWF_SegAccumE = 0;

			dAccumCS += nSegAccumCS / double(nGCEWF_AvgC);
			nSegAccumCS = 0;

		}
	}
}
//======== Wang-Landau ==========//

void Ising::SetWL(double * lnGe_p, unsigned long *geN_p, int size,
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

bool Ising::WLCheckBelowAVG() {
	double dAvg = uCurMove;
	dAvg /= nEn - 2;
	dAvg *= FlatRate;
	for (int n = 0; n < nEn; n++) {
		cout << GeN_p[n] << '\t';
	}
	cout << endl;

	//skip the unaccessible E
	if (GeN_p[0] < dAvg)
		return false;
	if (GeN_p[nEn - 1] < dAvg)
		return false;
	for (int n = 2; n < nEn - 2; n++) {
		if (GeN_p[n] < dAvg)
			return false;
	}
	return true;
}
void Ising::WLResetGeN() {
	uCurMove = 0;
	for (int n = 0; n < nEn; n++) {
		GeN_p[n] = 0;
		LnGe_p[n] -= LnGe_p[0];//save precisions as far as possible
	}

}
void Ising::WangLandauTrial(int Moves) {
	int CurStep = 0;
	int E1, E2;
	double dDeltLnGE_old2new;

	while (CurStep < Moves) {
		Delt();
		//cout << (2 * nN + nCurHami + nDeltaE) << endl;
		E1 = (2 * nN + nCurHami) / 4;
		E2 = (2 * nN + nCurHami + nDeltaE) / 4;

		dDeltLnGE_old2new = LnGe_p[E1] - LnGe_p[E2];
		if (dDeltLnGE_old2new < 0) {
			if (ran.Real() < exp(dDeltLnGE_old2new)) {
				LnGe_p[E2] += Lnf;
				GeN_p[E2] += 1;
				States_p[nCurTrial] *= -1;//< accepted
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
			States_p[nCurTrial] *= -1;//< accepted
			nCurHami += nDeltaE;
			//nAccepted = 1;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
}
//==========
bool Ising::WLCheckBelowAVGBound() {
	double dAvg = uCurMove;

	int lEi = (2 * nN + nLowEBound) / 4;
	int hEi = (2 * nN + nHighEBound) / 4;
	//assume (lEi >= 6)
	dAvg /= (nHighEBound - nLowEBound) / 4 + 1;
	dAvg *= FlatRate;
	for (int n = lEi; n <= hEi; n++) {
		if (GeN_p[n] < dAvg)
			return false;
	}
	return true;
}
void Ising::WLResetGeNBound() {
	uCurMove = 0;
	double reflnGe = LnGe_p[(2 * nN + nLowEBound) / 4];
	for (int n = 0; n < nEn; n++) {
		GeN_p[n] = 0;
		LnGe_p[n] -= reflnGe;//save precisions as far as possible
	}
}

void Ising::WangLandauTrialBound(int Moves) {
	int CurStep = 0;
	int E1, E2;
	double dDeltLnGE_old2new;

	while (CurStep < Moves) {
		Delt();//Assume nLowEBound<=nCurHami<=nHighEBound
		if ((nCurHami + nDeltaE) < nLowEBound || (nCurHami + nDeltaE)
				> nHighEBound)
			continue;

		E1 = (2 * nN + nCurHami) / 4;
		E2 = (2 * nN + nCurHami + nDeltaE) / 4;

		dDeltLnGE_old2new = LnGe_p[E1] - LnGe_p[E2];

		if (dDeltLnGE_old2new < 0) {
			if (ran.Real() < exp(dDeltLnGE_old2new)) {
				LnGe_p[E2] += Lnf;
				GeN_p[E2] += 1;
				States_p[nCurTrial] *= -1;//< accepted
				nCurHami += nDeltaE;
				nCurMag += nDeltaMag;
				//nAccepted = 0;
			} else {
				LnGe_p[E1] += Lnf;
				GeN_p[E1] += 1;
				//nAccepted = -1;
			}
		} else {
			LnGe_p[E2] += Lnf;
			GeN_p[E2] += 1;
			States_p[nCurTrial] *= -1;//< accepted
			nCurHami += nDeltaE;
			nCurMag += nDeltaMag;
			//nAccepted = 1;
		}
		uCurMove++;
		CurStep += 1;
	}
	Arg_io[AI_HAMI] = nCurHami;
	Arg_io[AI_MAG] = nCurMag;
}
void Ising::WarmupToBound() {
	InitConf(-1);//warm up
	Gauge();
	while (true) {
		Delt();
		States_p[nCurTrial] *= -1;//< accepted,T=0
		nCurHami += nDeltaE;
		//cout << nCurHami << endl;
		if ((nLowEBound <= nCurHami) && (nCurHami <= nHighEBound)) {
			break;
		}
	}
}
//========== HEAT Bath ===========//
void Ising::SetHeatBath(double T) {
	Arg_io[AI_TEMPER] = T;
	double h;
	int i;
	for (i = 0; i < 5; i++) //Table for heat bath
	{
		h = (-4 + 2 * i) / T;//for current spin=-1,J=1

		uLocalFieldPro[i] = (unsigned int) (TO32 * (exp(-h)
				/ (exp(h) + exp(-h))));//for j=1,map to uint32_t space,
#ifdef _DEBUG
		cout << h << '\t' << uLocalFieldPro[i] << endl;
#endif
	}
}
void Ising::HeatBathTrial(int Moves) {
	int LCField = 0, XTotalField = 0;
	unsigned int neib;

	for (int CurStep = 0; CurStep < Moves; CurStep++) {
		LCField = 0;
		nCurTrial = ran.Number(0, nN - 1);

		neib = nCurTrial + 1;//right
		if (neib % nDim1 == 0)
			neib -= nDim1;
		LCField += States_p[neib];

		neib = nCurTrial - 1;//left
		if (nCurTrial % nDim1 == 0)
			neib += nDim1;
		LCField += States_p[neib];

		neib = nCurTrial - nDim1;//up
		if (nCurTrial < nDim1)
			neib += nN;
		LCField += States_p[neib];

		neib = nCurTrial + nDim1;//down
		if (neib >= nN)
			neib -= nN;
		LCField += States_p[neib];

		//index in the probability table according to h=(4+local_field)/2
		if (ran.uNumber32() < uLocalFieldPro[(4 + LCField) >> 1]) {
			XTotalField += LCField * (States_p[nCurTrial] != -1);
			States_p[nCurTrial] = -1;
		} else {
			XTotalField -= LCField * (States_p[nCurTrial] != 1);
			States_p[nCurTrial] = 1;
		}

		uCurMove++;

	}
	nDeltaE = (XTotalField << 1);
	Arg_io[AI_DELT_HAMI] = nDeltaE;
	Arg_io[AI_HAMI] += Arg_io[AI_DELT_HAMI];
}

