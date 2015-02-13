/*
 ============================================================================
 Name        : RCEXPotts.cpp
 Description : Sampling in Potts model with rotational canonical Ensemble:beta(E)=beta0+alpha*(E-E0)/N.
 ============================================================================
 */

#include "Potts.h"
#include "nvector.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>


using namespace std;

class CRCEXPotts: public Potts {
protected:
	double *SE;
	double *E;
	double UN0, Beta0, Alpha;
	int nInitconf;
	double mean, std, mean1, std1;
	unsigned int uTau;//Time of AutoCorrelation (MCS/site)
	unsigned int uTeq, uTdrop;//Equilibrium at Teq
	unsigned int neff;
	double Beta;
public:
	char filename[120];
	char fileprefix[120];
	CRCEXPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		nInitconf = initconf;
		InitConf(nInitconf);
		E = new double[nEn];
		for (int n = 0; n < nEn; n++) {
			E[n] = nEn - n - 1;
			E[n] *= -1;
		}
		SE = new double[nEn];
	}

	void SetSE(double uN0, double beta0, double alpha) {
		UN0 = uN0;
		Beta0 = beta0;
		Alpha = alpha;
		double se0 = 0;
		for (int n = 0; n < nEn; n++) {
			if ((Beta0 + Alpha / nN * (E[n] - UN0 * nN)) < 0) {//avoid that beta is negative
				SE[n] = 0;
				se0 = E[n] * Beta0 + Alpha / nN * ((E[n] - UN0 * nN) * (E[n]
						- UN0 * nN)) / 2.0;
			} else {

				SE[n] = E[n] * Beta0 + Alpha / nN * ((E[n] - UN0 * nN) * (E[n]
						- UN0 * nN)) / 2.0 - se0;
			}
		}
		SetSefun(SE, nEn);
	}
	void RunSE_Ex(double uNt, double ExSigma, double CutSigma,
			unsigned long Samples, unsigned int Stride) {
		double duN = UN0 - uNt;
		double un0t;
		bool cont = false;//low(uNt)<-hight
		do {
			RunMetroplisSE(Samples, Stride);
			sprintf(filename, "%sQuanFull.txt", fileprefix);
			RunStatQuan(filename, Samples, Samples / 2);
			if (duN > 1.8 * std) {
				un0t=
				cont = true;
				if (uTdrop < Samples * 0.3) {
					Alpha -= Alpha * 0.05 * (Samples * 0.3 - uTdrop) / (Samples
							* 0.3);//go on
				} else if (uTdrop > Samples * 0.7) {
					Alpha += Alpha * 0.05 * (uTdrop - Samples * 0.7) / (Samples
							* 0.3);//redo
				}//else alpha+=0;
			} else {
				cont = false;
			}

		} while (cont);
	}
	void RunMetroplisSE(unsigned long Samples, unsigned int Stride) {
		sprintf(fileprefix, "rcepx%d_%d_%d~%d~%lf_%lf_%lf_%ld", nDim0, nDim1,
				int(cQ), nInitconf, UN0, Beta0, Alpha, Samples);
		sprintf(filename, "%sQuanFull.txt", fileprefix);
		std::ofstream ofile(filename);

		Gauge();
		Stride *= nN;

		//unsigned int QDis[cQ + 1];
		//double m1, m2;
		nvector<char> cons(States_p, nN);
		//unsigned long duS = Samples / 30;
		for (unsigned long s = 0; s < Samples; s++) {
			MetroplisTrialSe(Stride);
			/*m2 = QDistri(QDis, cQ + 1);
			 m1 = cQ * QDis[0];
			 m1 = (m1 / nN - 1) / (cQ - 1);
			 */
			ofile << setprecision(16) << Arg_io[AI_HAMI] << endl;
			/*if (s % duS == 0) {
			 sprintf(filename, "%sMidC%03ld_%lg_%lg_%lg.txt", fileprefix, s
			 / duS, Arg_io[AI_HAMI], m1, m2);
			 cons.dump_bin(filename);
			 }
			 */
		}
		ofile.close();

		sprintf(filename, "%sLastC.txt", fileprefix);//the last conf.
		cons.dump_bin(filename);

	}

	int RunStatQuan(char * file, unsigned long Samples, unsigned long rehead) {
		nvector<double> Esamples;
		Esamples.load(file, 0);
		if (Esamples.len < Samples) {
			cerr << "#Not so much data!" << endl;
			return -1;
		}

		unsigned int win = (unsigned int) (8 * pow(nN, 0.52));
		unsigned long uLags = Esamples.statis(500, win, uTeq, uTau);//~3*tau);
		uTdrop = uTeq + rehead;//Esamples.len / 10;//remove head 10%
		if (uTdrop > Samples) {
			cout << "uTdrop>Samples" << endl;
			uTdrop = uTeq;
		}

		if (uLags == 0)
			uTau = 1;

		neff = Esamples.mean_std(&mean, &std, false, uTdrop, uTau);
		Beta = Beta0 + Alpha * (mean / nN - UN0);
		Esamples.mean_std(&mean1, &std1, false, uTdrop, 1);

		cout << "#nInit\t" << setw(9) << "UN0"
				<< "\tB0\tAlpha\tWin\tN\tTeq\tnLags\tTau\tNeff\tmean\tstd\tBeta\tmean1\tstd1"
				<< endl << nInitconf << '\t' << setw(9) << setprecision(16)
				<< UN0 << '\t' << Beta0 << '\t' << Alpha << '\t' << win << '\t'
				<< Esamples.len << '\t' << uTdrop << '\t' << uLags << '\t'
				<< uTau << '\t' << neff << '\t' << mean << '\t' << std << '\t'
				<< Beta << '\t' << mean1 << '\t' << std1 << endl;
		return 0;
	}

	~ CRCEXPotts() {
		delete[] States_p;
		delete[] SE;
		delete[] E;
	}
};

int main(int argc, char *argv[]) {

	if (argc == 13) {
		clock_t __start = clock();
		CRCEXPotts
				p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
		p.SetSE(atof(argv[6]), atof(argv[7]), atof(argv[8]));
		p.RunSE_Ex(atof(argv[5]), atof(argv[11]), atof(argv[12]),
				atoi(argv[9]), atoi(argv[10]));
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.fileprefix << std::endl;
	} else if (argc == 3) {
		int L0, L1, Q, initconf;//!!!
		unsigned long Samples;
		double UN0, B0, alpha;
		int sv = sscanf(argv[1], "rcepx%d_%d_%d~%d~%lf_%lf_%lf_%ldQuan", &L0,
				&L1, &Q, &initconf, &UN0, &B0, &alpha, &Samples);
		if (sv != 8) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CRCEXPotts p(L0, L1, Q, initconf);
		p.SetSE(UN0, B0, alpha);//UN0,B0,alpha
		p.RunStatQuan(argv[1], Samples, atol(argv[2]));
	} else {
		cout
				<< "<RCEXPotts> l0,l1,Q,initconf,UN0,B0,alpha,samples,stride,dUN,Nrep"
				<< endl;//sampling at first time
		cout << "<RCEXPotts> rcepx%d_%d_%d~%d~%lf_%lf_%lf_%ldQuan*.txt nLags"
				<< endl;//statistics
		return -1;
	}
	return 0;
}

