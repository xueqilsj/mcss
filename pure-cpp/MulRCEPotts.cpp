/*
 ============================================================================
 Name        : MulRCEPotts.cpp
 Description : Multi-generalized canonical simulation :beta(E)=beta0+alpha*(E-E0)/N.
 ============================================================================
 */

#include "Potts.h"
#include "nvector.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

#define HIST_CUT 20

using namespace std;

double BetaFun(double *Es, double *Betas, double * Sigmas, double *Alphas,
		int Ns, double energy) {

	double beta, dSlope_t;
	int k = 0;
	if (energy <= Es[k]) {
		beta = Betas[k] + Alphas[k] * (energy - Es[k]);
		if (beta < 0)
			beta = 0;
	} else {
		while (k != Ns - 1 && energy > Es[k]) {
			k++;
		}
		if (energy > Es[Ns - 1]) {
			beta = Betas[Ns - 1] + Alphas[Ns - 1] * (energy - Es[Ns - 1]);
		} else {
			dSlope_t = (Betas[k + 1] - Betas[k]) / (Es[k + 1] - Es[k + 1]);
			beta = Betas[k] + dSlope_t * (energy - Betas[k]);
		}
	}
	return beta;
}
#define MAX_EX 40
class CMulRCEPotts: public Potts {
protected:
	double *SE, *E, *H, *WE;
	double Es[MAX_EX], Betas[MAX_EX], Sigmas[MAX_EX], Alphas[MAX_EX];
	int nInitconf, pn;
public:
	char fileprefix[120];
	char filename[120];
	CMulRCEPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		nInitconf = initconf;
		InitConf(nInitconf);
		E = new double[nEn];
		H = new double[nEn];
		for (int n = 0; n < nEn; n++) {
			E[n] = nEn - n - 1;
			E[n] *= -1;
		}
		SE = new double[nEn];
		WE = new double[nEn];
		pn = 0;
	}

	void PushB_E(double U0, double beta0, double alpha, double sigma) {
		Es[pn] = U0, Betas[pn] = beta0, Sigmas[pn] = sigma, Alphas[pn] = alpha;
		pn++;
		SE[0] = 0;
		WE[0] = 1;
		for (int n = 1; n < nEn; n++) {
			SE[n] = SE[n - 1] + BetaFun(Es, Betas, Sigmas, Alphas, pn, E[n]);
		}
		double se0 = 0;

		for (int n = 0; n < nEn; n++) {
			if ((beta0 + alpha / nN * (E[n] - U0)) < 0) {//avoid that beta is negative
				WE[n] = 1;
				se0 = E[n] * beta0 + alpha / nN * ((E[n] - U0) * (E[n] - U0))
						/ 2.0;
			} else {
				SE[n] = E[n] * beta0 + alpha / nN * ((E[n] - U0) * (E[n] - U0))
						/ 2.0 - se0;
				WE[n] = exp(-SE[n]);
			}
		}

		SetSefun(SE, nEn);
	}
	void RunMetroplisSENoConf(double UN0, double Beta0, double Alpha,
			double sigma) {

		sprintf(fileprefix, "mulrcep%d_%d_%d~%d~%lf_%lf_%lf", nDim0, nDim1,
				int(cQ), nInitconf, UN0, Beta0, Alpha);
		sprintf(filename, "%sStat.txt", fileprefix);
		std::ofstream ofile(filename);

		Gauge();
		double mean, mean2;
		unsigned long Samples;
		mean = UN0 * nN;
		mean2 = 0;

		cout << mean << '\t' << sigma << endl;
		PushB_E(mean, Beta0, Alpha, sigma);
		for (int n = 0; n < nEn; n++) {
			ofile << setprecision(16) << SE[n] << '\t' << H[n] << endl;
		}
		ofile << endl << endl;

		for (int n = 0; n < nEn; n++) {
			H[n] = 0;
		}

		Samples = 10000;
		MetroplisTrialSe(Samples / 10);

		int nk;
		for (unsigned long s = 0; s < Samples; s++) {
			MetroplisTrialSe(nN);
			nk = nEn - 1 + nCurHami;
			H[nk] += 1;
		}

		mean += E[nk];
		mean2 += E[nk] * E[nk];

		mean /= Samples;
		sigma = sqrt(mean2 / Samples - mean * mean);

		ofile.close();

	}

	int RunStatQuan(char * file, unsigned long Samples, unsigned long rehead) {
		return 1;
	}

	~ CMulRCEPotts() {
		delete[] States_p;
		delete[] SE;
		delete[] WE;
		delete[] E;
		delete[] H;
	}
};

int main(int argc, char *argv[]) {
	if (argc == 11) {
		clock_t __start = clock();
		CMulRCEPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(
				argv[4]));
		p.RunMetroplisSENoConf(atof(argv[5]), atof(argv[6]), atof(argv[7]),
				atoi(argv[8]));

		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.fileprefix << std::endl;

	} else {
		cout << "<MulRCEPotts> l0,l1,Q,initconf,UN0,B0,alpha" << endl;//sampling at first time
		return -1;
	}
	return 0;
}

