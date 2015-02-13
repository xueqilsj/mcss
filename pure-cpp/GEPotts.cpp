/*
 ============================================================================
 Name        : GEPotts.cpp
 Description : Sampling in Potts model with Generalized Ensemble:A part of LnGe (i.e. S(E)) plus two artificial boundary beta(E).
 ============================================================================
 */

#include "Potts.h"
#include "nvector.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;

class CGEPotts: public Potts {
	double *SE;
	double *E;
public:
	char filename[128];
	CGEPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(initconf);
		E = new double[nEn];
		for (int n = 0; n < nEn; n++) {
			E[n] = nEn - n - 1;
			E[n] *= -1;
		}
		SE = new double[nEn];
	}
	bool SetBetaE(char *file, double lt, double ht, int startN) {

		nvector<double> fBE;
		if (fBE.load(file, 0) < 0)
			return false;
		//cout << fBE << endl;

		int nMax = fBE.len;
		if (nMax + startN > nEn)
			nMax = nEn - startN;

		double beta;
		if (lt < 0) {
			if (nMax < 2)
				return false;
			beta = fBE.data[0];
		} else {
			beta = 1.0 / lt;
		}
		SE[0] = 0;
		for (int n = 1; n < startN; n++) {
			SE[n] = SE[n - 1] + beta;
			//cout << n - 1 << ":\t" << SE[n - 1] << endl;
		}

		for (int n = startN; n < startN + nMax; n++) {
			SE[n] = SE[n - 1] + fBE.data[n - startN];
			//cout << n << "=\t" << SE[n] << endl;
		}

		if (ht < 0) {
			if (nMax < 2)
				return false;
			beta = fBE.data[nMax - 1];
		} else {
			beta = 1.0 / ht;
		}

		for (int n = startN + nMax; n < nEn; n++) {
			SE[n] = SE[n - 1] + beta;
			//cout << n << "v\t" << SE[n] << endl;
		}
		SetSefun(SE, nEn);
		for (int n = 0; n < nEn; n++)
			cout << setprecision(15) << SE[n] << endl;
		return true;
	}
	bool SetSE(char *file, double lt, double ht, int startN) {
		//lt:low temperature,ht
		//negative :extend
		nvector<double> fSE;
		if (fSE.load(file, 0) < 0)
			return false;
		//cout << fSE << endl;

		int nMax = fSE.len;
		if (nMax + startN > nEn)
			nMax = nEn - startN;

		for (int n = startN; n < startN + nMax; n++) {
			SE[n] = fSE.data[n - startN];
			//cout << n << "=\t" << SE[n] << endl;
		}

		double beta;
		if (lt < 0) {
			if (nMax < 2)
				return false;
			beta = fSE.data[1] - fSE.data[0];
		} else {
			beta = 1.0 / lt;
		}

		for (int n = startN; n > 0; n--) {
			SE[n - 1] = SE[n] - beta;
			//cout << n - 1 << ":\t" << SE[n - 1] << endl;
		}

		if (ht < 0) {
			if (nMax < 2)
				return false;
			beta = fSE.data[nMax - 1] - fSE.data[nMax - 2];
		} else {
			beta = 1.0 / ht;
		}

		for (int n = startN + nMax; n < nEn; n++) {
			SE[n] = SE[n - 1] + beta;
			//cout << n << "v\t" << SE[n] << endl;
		}
		SetSefun(SE, nEn);
		for (int n = 0; n < nEn; n++)
			cout << setprecision(15) << SE[n] << endl;
		return true;
	}
	int RunMetroplisTrialSe(char *file, double lt, double ht, int startN,
			unsigned long Samples, unsigned int Stride, int nConfs) {

		if (strstr(file, "SE.txt") != NULL) {
			if (!SetSE(file, lt, ht, startN)) {
				cerr << "SetSE() Error!" << endl;
				return -1;
			}
		} else if (strstr(file, "BE.txt") != NULL) {
			if (!SetBetaE(file, lt, ht, startN)) {
				cerr << "SetBE() Error!" << endl;
				return -1;
			}
		} else {
			return -1;
		}

		long nl = strlen(file);
		strncpy(filename, file, nl - 4);//extern char *strncpy(char *dest, char *src, int n);

		sprintf(filename + nl - 4, "_%d~%g_%g_%ld_%dTraj.txt", nInitconf, lt,
				ht, Samples, Stride);
		Gauge();

		std::ofstream ofile(filename);
		Stride *= nN;
		nvector<double> TR(Samples);
		nvector<double> MTR(Samples);
		unsigned int *QDis = new unsigned int[cQ + 1];
		double m1, m2;
		unsigned long snaps;
		if (nConfs <= 0)
			snaps = Samples + 1;
		else
			snaps = Samples / nConfs;

		for (unsigned long CurSample = 0; CurSample < Samples; CurSample++) {
			MetroplisTrialSe(Stride);
			m2 = QDistri(QDis, cQ + 1);
			m1 = cQ * QDis[0];
			m1 = (m1 / nN - 1) / (cQ - 1);
			ofile << setprecision(16) << Arg_io[AI_HAMI] << '\t'
					<< setprecision(6) << m1 << '\t' << m2 << endl;
			if (CurSample % snaps == 0) {
				sprintf(filename + nl - 4, "_%d~%g_%g_%ld_%d~%03ld_%d.conf",
						nInitconf, lt, ht, Samples, Stride / nN, CurSample
								/ snaps, nCurHami);
				DumpConf(filename);
			}
		}
		delete[] QDis;
		ofile.close();
		return 0;
	}
	~ CGEPotts() {
		delete[] States_p;
		delete[] SE;
		delete[] E;
	}
};
int main(int argc, char *argv[]) {
	if (argc == 8) {
		clock_t __start = clock();
		int L0, L1, Q, startN;//!!!
		int sv = sscanf(argv[1], "gep%d_%d_%d~%d", &L0, &L1, &Q, &startN);
		if (sv != 4) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CGEPotts p(L0, L1, Q, atoi(argv[2]));

		p.RunMetroplisTrialSe(argv[1], atof(argv[3]), atof(argv[4]), startN,
				atoi(argv[5]), atoi(argv[6]), atoi(argv[7]));//SE.txt

		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;
	} else {
		cout
				<< "<GEPotts> gep%d_%d_%d~%d*SE.txt|gep%d_%d_%d~%d*BE.txt initf lowT highT samples stride nConf"
				<< endl;//L1,L2,Q,StartAtEnergyArrayIndex,lowT<0 or highT<0 straight line with slop same as boundary.
		return -1;
	}
	return 0;
}

