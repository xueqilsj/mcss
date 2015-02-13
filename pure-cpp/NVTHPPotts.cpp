/*
 ============================================================================
 Name        : HPNVTPotts.cpp
 Description : Sampling in NVT Potts model with Heat Path Algorithm
 ============================================================================
 */

#include "Potts.h"
#include "nvector.h"

#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace std;

class CNVTHPPotts: public Potts {
	double *E;
	int nInitconf;
public:
	char filename[80];
	CNVTHPPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		nInitconf = initconf;
		InitConf(nInitconf);
		E = new double[nEn];
		for ( int n = 0; n < nEn; n++) {
			E[n] = nEn - n - 1;
			E[n] *= -1;
		}
		QDeltaEi = new int[Q];

	}
	void RunMetroplis(double temp, unsigned long Samples, unsigned int Stride) {

		SetHeatPath(temp, QDeltaEi, cQ);
		for (int mE = 0; mE <= 4; mE++) {//normalization
			cout << _dDeltaExp[mE] << endl;
		}
		Gauge();

		sprintf(filename, "hpnvtp%d_%d_%d~%d~%lg_%ld_%dQuanFull.txt", nDim0,
				nDim1, int(cQ), nInitconf, temp, Samples, Stride);

		Stride *= nN;
		nvector<double> TR(Samples);
		nvector<double> MTR(Samples);

		unsigned long CurSample = 0, otSamples = 0;
		do {
			if (Samples - CurSample > 4000000)//4e6~2.2hours
				otSamples += 4000000;
			else
				otSamples = Samples;

			for (; CurSample < otSamples; CurSample++) {
				HeatPathTrial(1);

				TR.data[CurSample] = Arg_io[AI_HAMI];
				MTR.data[CurSample] = Arg_io[AI_MAG];
			}
			nv_dump(TR, MTR, filename);//save sample per 1e6
		} while (CurSample < Samples);

	}
	int RunStatQuan(char * file, unsigned long Samples, unsigned long nLags) {
		nvector<double> Esamples;
		Esamples.load(file, 0);
		if (Esamples.len < Samples) {
			cerr << "#Not so much data!" << endl;
			return -1;
		}

		unsigned int uTau;//Time of AutoCorrelation (MCS/site)
		unsigned int uTeq;//Equilibrium at Teq
		unsigned int win = (unsigned int) (8 * pow(nN, 0.52));
		unsigned long uLags = Esamples.statis(nLags, win, uTeq, uTau);//~3*tau);
		if (uLags == 0)
			uTau = 1;
		double mean, std, mean1, std1;
		unsigned int neff = Esamples.mean_std(&mean, &std, false, uTeq, uTau);
		Esamples.mean_std(&mean1, &std1, false, uTeq, 1);

		cout
				<< "#nInit\tTemp\tWin\tN\tTeq\tnLags\tTau\tNeff\tmean\tstd\tmean1\tstd1"
				<< endl << nInitconf << '\t' << Arg_io[AI_TEMPER] << '\t'
				<< win << '\t' << Esamples.len << '\t' << uTeq << '\t' << uLags
				<< '\t' << uTau << '\t' << neff << '\t' << mean << '\t' << std
				<< '\t' << mean1 << '\t' << std1 << endl;
		return 0;
	}
	~ CNVTHPPotts() {
		delete[] States_p;
		delete[] E;
		delete[] QDeltaEi;
	}
};
int main(int argc, char *argv[]) {
	if (argc == 8) {
		clock_t __start = clock();
		CNVTHPPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
				atoi(argv[4]));
		p.RunMetroplis(atof(argv[5]), atoi(argv[6]), atoi(argv[7]));
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;
	} else if (argc == 3) {
		int L0, L1, Q, initconf;//!!!
		unsigned long Samples;
		unsigned int stride;
		double temp;
		int sv = sscanf(argv[1], "nvtp%d_%d_%d~%d~%lg_%ld_%dQuan", &L0, &L1,
				&Q, &initconf, &temp, &Samples, &stride);
		if (sv != 7) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CNVTHPPotts p(L0, L1, Q, initconf);
		p.Arg_io[p.AI_TEMPER] = temp;
		p.RunStatQuan(argv[1], Samples, atol(argv[2]));

	} else {
		cout << "<NVTHPPotts> dim0 dim0 Q initConf temperature samples stride"
				<< endl;//stride =MCSweep/site
		cout << "<NVTHPPotts> nvtp%d_%d_%d~%d~%lg_%ld_%dQuan*.txt nLags"
				<< endl;//nLags for AutoCorr
		return -1;
	}
	return 0;
}

