/*
 ============================================================================
 Name        : GCEPotts.cpp
 Description : Sampling in Potts model with rotational canonical Ensemble:beta(E)=beta0+alpha*(E-E0)/N.
 ============================================================================
 */

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

#include "Potts.h"
#include "nvector.h"

using namespace std;

class CGCEPotts: public Potts {
protected:
public:
	char filename[128];
	CGCEPotts() {
	}
	CGCEPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(initconf);
	}

	void RunMetroplisTrialGCE(unsigned long Samples, unsigned long Stride,
			long nConfs) {
		sprintf(filename,
				"gcep%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ldTraj.txt", nDim0,
				nDim1, int(cQ), nInitconf, dBeta, dAlpha, dUn, Samples, Stride);
		std::ofstream ofile(filename);
		ofile << setprecision(16);

		Gauge();
		Stride *= nN;

		unsigned long snaps;
		if (nConfs <= 0)
			snaps = Samples + 1;
		else
			snaps = Samples / nConfs;

		for (unsigned long CurSample = 1; CurSample <= Samples; CurSample++) {
			MetroplisTrialGCE(Stride);
			ofile << nCurHami << endl;//<< '\t' << m2 << endl;
			if (CurSample % snaps == 0) {
				sprintf(
						filename,
						"gcep%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ld~%03ld_%d.conf",
						nDim0, nDim1, int(cQ), nInitconf, dBeta, dAlpha, dUn,
						Samples, Stride / nN, CurSample / snaps, nCurHami);
				DumpConf(filename);
			}
		}
		ofile.close();

	}

	int RunStatQuan(char * file, unsigned long Samples, unsigned long rehead) {
		nvector<double> Esamples;
		Esamples.load(file, 0);
		if (Esamples.len < Samples) {
			cerr << "#Not so much data!" << endl;
			return -1;
		}
		unsigned int uTau;//Time of AutoCorrelation (MCS/site)
		unsigned int uTeq, uTdrop;//Equilibrium at Teq
		unsigned int win = Samples / 50;//=2% (unsigned int) (8 * pow(nN, 0.52));
		unsigned long uLags = Esamples.statis(500, win, uTeq, uTau);//~3*tau);
		uTdrop = rehead;//Esamples.len / 10;//remove head 10%
		if (uTdrop < uTeq)
			uTdrop = uTeq;

		if (uLags == 0)
			uTau = 1;

		double mean, std, mean1, std1, meanErr;
		Esamples.mean_std_err(&mean, &std, &meanErr, uTdrop, uTau);
		double Beta = dBeta + dAlpha * (mean / nN - dUn);
		unsigned long neff = Esamples.mean_std(&mean1, &std1, false, uTdrop,
				uTau);

		cout << "#nInit\t" << setw(9) << "dUn"
				<< "\tB0\tAlpha\tWin\tN\tTeq\tnLags\tTau\tNeff\tmean\tstd\tBeta\tMerr\tmean1\tstd1"
				<< endl << nInitconf << '\t' << setw(9) << setprecision(16)
				<< dUn << '\t' << dBeta << '\t' << dAlpha << '\t' << win
				<< '\t' << Esamples.len << '\t' << uTdrop << '\t' << uLags
				<< '\t' << uTau << '\t' << neff << '\t' << mean << '\t' << std
				<< '\t' << Beta << '\t' << '\t' << meanErr << '\t' << mean1
				<< '\t' << std1 << endl;
		return 0;
	}

	~ CGCEPotts() {
		delete[] States_p;
	}
};

int main(int argc, char *argv[]) {
	int L0, L1, Q, initconf;//!!!
	unsigned long Samples;
	unsigned long Stride;
	double un, B0, alpha;
	clock_t __start = clock();

	if (argc == 11) {
		CGCEPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
		p.dBeta = atof(argv[5]);
		p.dAlpha = atof(argv[6]);
		p.dUn = atof(argv[7]);
		p.RunMetroplisTrialGCE(atol(argv[8]), atol(argv[9]), atol(argv[10]));
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;
	} else if (argc == 8) {
		CGCEPotts p;
		if (p.FromConf(argv[1])) {
			p.dBeta = atof(argv[2]);
			p.dAlpha = atof(argv[3]);
			p.dUn = atof(argv[4]);
			p.RunMetroplisTrialGCE(atol(argv[5]), atol(argv[6]), atol(argv[7]));
			cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
					<< p.filename << std::endl;
		}
	} else if (argc == 3) {
		int sv = sscanf(argv[1], "gcep%d_%d_%d~%d~%lg_%lg_%lg_%ld_%ldTraj",
				&L0, &L1, &Q, &initconf, &B0, &alpha, &un, &Samples, &Stride);
		if (sv != 9) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CGCEPotts p(L0, L1, Q, initconf);
		p.dUn = un;
		p.dBeta = B0;
		p.dAlpha = alpha;
		p.RunStatQuan(argv[1], Samples, atol(argv[2]));
	} else {
		cout << "<GCEPotts> l0 l1 Q initconf B0 alpha dUn samples stride nConf"
				<< endl;//sampling at first time
		cout
				<< "<GCEPotts> *%d_%d_%d~*Conf.txt B0 alpha dUn samples stride nConf"
				<< endl;//statistics
		cout
				<< "<GCEPotts> gcep%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ldTraj*.txt nLags"
				<< endl;//statistics
		return -1;
	}
	return 0;
}

