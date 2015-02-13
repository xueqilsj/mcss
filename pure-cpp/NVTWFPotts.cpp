/*
 ============================================================================
 Name        : WFPotts.cpp
 Description : Sampling in Potts model with Wollf Cluster method
 ============================================================================
 */
#include "Potts.h"
#include "nvector.h"

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <time.h>

using namespace std;

class CWFPotts: public Potts {

public:
	char filename[128];
	CWFPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(initconf);
	}
	void RunMetroplis(double temp, unsigned long Samples, unsigned int WFsteps,
			int nConfs) {
		if (temp <= 0) {
			dPadd = 0.5;
			Arg_io[AI_TEMPER] = log(2.0);
		} else {
			dPadd = 1 - exp(-temp);
			Arg_io[AI_TEMPER] = temp;
		}

		int *stack = new int[nN];
		ClusterStack_p = stack;
		Gauge();

		sprintf(filename, "wfp%d_%d_%d~%d~%.15lg_%ld_%dTraj.txt", nDim0, nDim1,
				int(cQ), nInitconf, temp, Samples, WFsteps);
		std::ofstream ofile(filename);

		unsigned long snaps;
		if (nConfs <= 0)
			snaps = Samples + 1;
		else
			snaps = Samples / nConfs;

		ofile << setprecision(16);
		for (unsigned long CurSample = 1; CurSample <= Samples; CurSample++) {
			WolffCluster(WFsteps);
			//hs = Arg_io[AI_HAMI];
			Gauge();//slow down!!
			if (CurSample % snaps == 0) {
				sprintf(filename, "wfp%d_%d_%d~%d~%.15lg_%ld_%d~%03ld_%d.conf",
						nDim0, nDim1, int(cQ), nInitconf, temp, Samples,
						WFsteps, CurSample / snaps, nCurHami);
				DumpConf(filename);
			}
			ofile << nCurHami << endl;
		}
		ofile.close();
		delete[] stack;
	}
	int RunStatQuan(char * file, unsigned long Samples, unsigned long nLags) {
		dPadd = 1 - exp(-Arg_io[AI_TEMPER]);
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
				<< "#nInit\tTemp\tPadd\tWin\tN\tTeq\tnLags\tTau\tNeff\tmean\tstd\tmean1\tstd1"
				<< endl << nInitconf << '\t' << Arg_io[AI_TEMPER] << '\t'
				<< dPadd << '\t' << win << '\t' << Esamples.len << '\t' << uTeq
				<< '\t' << uLags << '\t' << uTau << '\t' << neff << '\t'
				<< mean << '\t' << std << '\t' << mean1 << '\t' << std1 << endl;
		return 0;
	}
	~ CWFPotts() {
		delete[] States_p;

	}
};
int main(int argc, char *argv[]) {
	if (argc == 9) {
		clock_t __start = clock();
		CWFPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
		p.RunMetroplis(atof(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(
				argv[8]));
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;
	} else if (argc == 3) {
		int L0, L1, Q, initconf;//!!!
		unsigned long Samples;
		unsigned int stride;
		double temp;
		int sv = sscanf(argv[1], "wfp%d_%d_%d~%d~%lg_%ld_%dTraj", &L0, &L1, &Q,
				&initconf, &temp, &Samples, &stride);
		if (sv != 7) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CWFPotts p(L0, L1, Q, initconf);

		if (temp <= 0) {
			p.Arg_io[p.AI_TEMPER] = 1.0 / log(2.0);
		} else {
			p.Arg_io[p.AI_TEMPER] = temp;
		}

		p.RunStatQuan(argv[1], Samples, atol(argv[2]));

	} else {
		cout << "<NVTWFPotts> dim0 dim0 Q initConf beta samples WFsteps nConf"
				<< endl;//WFstep :per Cluster
		cout << "<NVTWFPotts> wfp%d_%d_%d~%d~%lg_%ld_%dTraj*.txt nLags" << endl;//nLags for AutoCorr
		return -1;
	}
	return 0;

}

