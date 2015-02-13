/*
 ============================================================================
 Name        : NVTPotts.cpp
 Description : Sampling in NVT Potts model with Metroplis Algorithm
 ============================================================================
 */

#include "Potts.h"
#include "nvector.h"

#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace std;

class CNVTPotts: public Potts {
	double *E;
public:
	char filename[128];
	CNVTPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(initconf);
	}
	void RunMetroplisTrialBT(double temp, unsigned long Samples,
			unsigned int Stride, int nConfs) {
		//SetTemperature(temp);
		Arg_io[AI_TEMPER] = temp;
		Gauge();

		sprintf(filename, "nvtp%d_%d_%d~%d~%lf_%ld_%dTraj.txt", nDim0, nDim1,
				int(cQ), nInitconf, temp, Samples, Stride);
		Stride *= nN;
		std::ofstream ofile(filename);
		ofile << setprecision(16);
		unsigned long snaps;
		if (nConfs <= 0)
			snaps = Samples + 1;
		else
			snaps = Samples / nConfs;

		for (unsigned long CurSample = 1; CurSample <= Samples; CurSample++) {
			MetroplisTrialBT(Stride);
			ofile << Arg_io[AI_HAMI] << endl;//<< '\t' << m2 << endl;
			if (CurSample % snaps == 0) {
				sprintf(filename, "nvtp%d_%d_%d~%d~%lf_%ld_%d~%03ld_%d.conf",
						nDim0, nDim1, int(cQ), nInitconf, temp, Samples, Stride
								/ nN, CurSample / snaps, nCurHami);
				DumpConf(filename);
			}
		}
		ofile.close();
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
	~ CNVTPotts() {
		delete[] States_p;
	}
};
int main(int argc, char *argv[]) {
	if (argc == 9) {
		clock_t __start = clock();
		CNVTPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
		p.RunMetroplisTrialBT(atof(argv[5]), atoi(argv[6]), atoi(argv[7]),
				atoi(argv[8]));
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;
	} else if (argc == 3) {
		int L0, L1, Q, initconf;//!!!
		unsigned long Samples;
		unsigned int stride;
		double temp;
		int sv = sscanf(argv[1], "nvtp%d_%d_%d~%d~%lf_%ld_%dTraj", &L0, &L1,
				&Q, &initconf, &temp, &Samples, &stride);
		if (sv != 7) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CNVTPotts p(L0, L1, Q, initconf);
		p.Arg_io[p.AI_TEMPER] = temp;
		p.RunStatQuan(argv[1], Samples, atol(argv[2]));

	} else {
		cout << "<NVTPotts> dim0 dim0 Q initConf beta samples stride nConf"
				<< endl;//stride =MCSweep/site
		cout << "<NVTPotts> nvtp%d_%d_%d~%d~%lf_%ld_%dTraj.txt nLags" << endl;//nLags for AutoCorr
		return -1;
	}
	return 0;
}

