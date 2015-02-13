/*
 ============================================================================
 Name        : RCEPotts.cpp
 Description : Sampling in Potts model with rotational canonical Ensemble:beta(E)=beta0+alpha*(E-E0)/N.
 ============================================================================
 */

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

#include "Potts.h"
#include "nvector.h"

#define SAMPLES_ONCE 20

using namespace std;

class CRCEPotts: public Potts {
protected:
	double *SE;
	double *E;
public:
	char filename[128];
	CRCEPotts(int dim0, int dim1, char Q, int initconf = 0) {
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
	~ CRCEPotts() {
		delete[] States_p;
		delete[] SE;
		delete[] E;
	}

	void SetSE(double beta0, double alpha, double uN0) {
		dUn = uN0;
		dBeta = beta0;
		dAlpha = alpha;
		double se0 = 0;
		for (int n = 0; n < nEn; n++) {
			if ((dBeta + dAlpha / nN * (E[n] - dUn * nN)) < 0) {//avoid that beta is negative
				SE[n] = 0;
				se0 = E[n] * dBeta + dAlpha / nN * ((E[n] - dUn * nN) * (E[n]
						- dUn * nN)) / 2.0;
			} else {

				SE[n] = E[n] * dBeta + dAlpha / nN * ((E[n] - dUn * nN) * (E[n]
						- dUn * nN)) / 2.0 - se0;
			}
		}
		SetSefun(SE, nEn);
	}
	void RunMetroplisTrialSe(unsigned long Samples, unsigned int Stride,
			int nConfs) {
		sprintf(filename, "rcep%d_%d_%d~%d~%lf_%lf_%.14lf_%ld_%dTraj.txt",
				nDim0, nDim1, int(cQ), nInitconf, dBeta, dAlpha, dUn, Samples,
				Stride);
		std::ofstream ofile(filename);

		Gauge();
		Stride *= nN;

		unsigned long snaps;
		if (nConfs <= 0)
			snaps = Samples + 1;
		else
			snaps = Samples / nConfs;

		ofile << setprecision(16);
		for (unsigned long CurSample = 1; CurSample <= Samples; CurSample++) {
			MetroplisTrialSe(Stride);
			ofile << nCurHami << endl;//<< '\t' << m2 << endl;
			if (CurSample % snaps == 0) {
				sprintf(
						filename,
						"rcep%d_%d_%d~%d~%lf_%lf_%.14lf_%ld_%d~%03ld_%d.conf",
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
		unsigned int win = (unsigned int) (8 * pow(nN, 0.52));
		unsigned long uLags = Esamples.statis(500, win, uTeq, uTau);//~3*tau);
		uTdrop = rehead;//Esamples.len / 10;//remove head 10%
		if (uTdrop < uTeq)
			uTdrop = uTeq;

		if (uLags == 0)
			uTau = 1;
		double mean, std, mean1, std1;

		unsigned int neff = Esamples.mean_std(&mean, &std, false, uTdrop, uTau);
		double Beta = dBeta + dAlpha * (mean / nN - dUn);
		Esamples.mean_std(&mean1, &std1, false, uTdrop, 1);

		cout << "#nInit\t" << setw(9) << "dUn"
				<< "\tB0\tAlpha\tWin\tN\tTeq\tnLags\tTau\tNeff\tmean\tstd\tBeta\tmean1\tstd1"
				<< endl << nInitconf << '\t' << setw(9) << setprecision(16)
				<< dUn << '\t' << dBeta << '\t' << dAlpha << '\t' << win
				<< '\t' << Esamples.len << '\t' << uTdrop << '\t' << uLags
				<< '\t' << uTau << '\t' << neff << '\t' << mean << '\t' << std
				<< '\t' << Beta << '\t' << mean1 << '\t' << std1 << endl;
		return 0;
	}

};

int main(int argc, char *argv[]) {

	if (argc == 11) {
		clock_t __start = clock();
		CRCEPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
		p.SetSE(atof(argv[5]), atof(argv[6]), atof(argv[7]));
		p.RunMetroplisTrialSe(atoi(argv[8]), atoi(argv[9]), atoi(argv[10]));

		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;
	} else if (argc == 5) {
		clock_t __start = clock();
		int L0, L1, Q, initconf;//!!!
		unsigned long Samples;
		double dUn, B0, alpha;
		if (sscanf(argv[1], "rcep%d_%d_%d~%d~%lf_%lf_%lf_%ld", &L0, &L1, &Q,
				&initconf, &B0, &alpha, &dUn, &Samples) != 8) {
			cout << "Error input:" << argv[1] << endl;
			return -1;
		} else if (strstr(argv[1], "LastC.txt") != NULL) {
			//"<RCEPotts> rcep%d_%d_%d~%d~%lf_%lf_%lf_%ldlastC*.txt samples,stride,saveConf"
			CRCEPotts p(L0, L1, Q, initconf);
			nvector<char> cons(p.States_p, p.nN);
			cons.load_bin(argv[1], p.nN, 0, false);//set last conf

			p.SetSE(B0, alpha, dUn);//dUn,B0,alpha
			p.RunMetroplisTrialSe(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));

			cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
					<< p.filename << std::endl;
		} else {
			cout << "Unknown input:" << argv[1] << endl;
			return -1;
		}

	} else if (argc == 3) {
		int L0, L1, Q, initconf;//!!!
		unsigned long Samples;
		double dUn, B0, alpha;
		int sv = sscanf(argv[1], "rcep%d_%d_%d~%d~%lf_%lf_%lf_%ldTraj", &L0,
				&L1, &Q, &initconf, &B0, &alpha, &dUn, &Samples);
		if (sv != 8) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CRCEPotts p(L0, L1, Q, initconf);
		p.SetSE(B0, alpha, dUn);//dUn,B0,alpha
		p.RunStatQuan(argv[1], Samples, atol(argv[2]));
	} else {
		cout << "<RCEPotts> l0 l1 Q initconf B0 alpha dUn samples stride nConf"
				<< endl;//sampling at first time
		cout
				<< "<RCEPotts> rcep%d_%d_%d~%d~%lf_%lf_%lf_%ldlastC*.txt samples stride nConf"
				<< endl;//start from last time

		cout << "<RCEPotts> rcep%d_%d_%d~%d~%lf_%lf_%lf_%ldTraj*.txt nLags"
				<< endl;//statistics
		return -1;
	}
	return 0;
}

