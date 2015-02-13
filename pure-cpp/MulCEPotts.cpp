/*
 ============================================================================
 Name        : MulCEPotts.cpp
 Description : Multi(generalized)canonical Ensemble simulation:
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

class CMulCEPotts: public Potts {
protected:
	double *SE;
	double *E;
	double *H;
public:
	char fileprefix[120];
	char filename[120];
	CMulCEPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(initconf);
		E = new double[nEn];
		H = new double[nEn];
		for (int n = 0; n < nEn; n++) {
			E[n] = nEn - n - 1;
			E[n] *= -1;
			H[n] = 0;
		}
		SE = new double[nEn];
	}
	void SmoothH(int minEn, int maxEn, int win) {
		double ek;
		for (int k = minEn + win; k < maxEn - win; k++) {
			ek = 0;
			for (int j = -win; j <= win; j++)
				ek += H[k + j];
			H[k] = ek / (2 * win + 1);
		}
	}
	void RunMetroplisSE_MulRCE2(double un0, double beta0, double alpha,
			double ENend, int iter) {
		double U0 = un0 * nN;
		Arg_io[AI_TEMPER] = beta0;
		for (int n = 0; n < nEn; n++) {
			SE[n] = E[n] * Arg_io[AI_TEMPER] + alpha / nN * ((E[n] - U0)
					* (E[n] - U0)) / 2.0;
		}
		SetSefun(SE, nEn);

		sprintf(fileprefix, "mulcep%d_%d_%d~%d~%lf_%lf_%lf_%lf_%d", nDim0,
				nDim1, int(cQ), nInitconf, un0, Arg_io[AI_TEMPER], alpha,
				ENend, iter);
		sprintf(filename, "%sStat.txt", fileprefix);
		std::ofstream ofile(filename);

		Gauge();
		unsigned long Samples = 10000;

		MetroplisTrialSe(Samples / 10);
		for (unsigned long s = 0; s < Samples; s++) {
			MetroplisTrialSe(nN);
			H[nEn - 1 + nCurHami] += 1;
		}
		for (int n = 0; n < nEn; n++) {
			ofile << setprecision(16) << SE[n] << '\t' << H[n] << endl;
		}
		ofile << endl << endl;

		double Eend = long(ENend * nN);
		double mean = 0, mean2 = 0;
		double mBeta, b0;
		int EminN = -1, EmaxN, EmaxN2;

		for (int n = 0; n < nEn; n++) {
			mean += H[n] * E[n];
			mean2 += H[n] * E[n] * E[n];
		}
		mean /= Samples;
		mean2 /= Samples;
		double sig = sqrt(mean2 - mean * mean);

		EmaxN = nEn - 1 + int(mean);
		EminN = nEn - 1 + int(mean - 1 * sig);

		do {
			if (Eend >= E[EminN]) {
				if (H[EmaxN] >= HIST_CUT)
					iter--;
				EminN = nEn - 1 + int(Eend);
			}
			for (int n = EmaxN; n >= EminN; n--) {
				if (H[n] >= HIST_CUT) {
					EmaxN2 = n;
					break;
				}
			}

			b0 = (H[EminN] > 0) ? H[EminN] : HIST_CUT;
			for (int n = EminN; n <= EmaxN2; n++) {
				if (H[n] >= b0) {
					SE[n] += log(H[n] / H[EmaxN2]);
				}
			}

			int k = 1;
			while (H[EminN + k] < HIST_CUT && k < EmaxN - EminN)
				k++;
			mBeta = (SE[EminN + k] - SE[EminN]) / k;

			for (int n = EminN; n > 0; n--) {
				b0 = mBeta + alpha / nN * (E[n] - E[EminN]);
				if (b0 < 0)
					b0 = 0;
				SE[n - 1] = SE[n] - b0;
			}

			for (int n = 0; n < nEn; n++) {
				H[n] = 0;
			}

			Samples = (EmaxN - EminN) * HIST_CUT * nDim1;
			cout << E[EminN] << '\t' << E[EmaxN2] << '\t' << E[EmaxN] << '\t'
					<< sig << '\t' << H[EminN] << '\t' << Samples << '\t'
					<< iter << endl;

			MetroplisTrialSe(Samples / 10);
			for (unsigned long s = 0; s < Samples; s++) {
				MetroplisTrialSe(nN);
				H[nEn - 1 + nCurHami] += 1;
			}

			//SmoothH(EminN, EmaxN, 4);

			for (int n = 0; n <= EminN; n++) {
				mean += H[n] * E[n];
				mean2 += H[n] * E[n] * E[n];
			}
			mean /= Samples;
			mean2 /= Samples;
			sig = sqrt(mean2 - mean * mean);
			EminN = nEn - 1 + int(mean - 1 * sig);

			for (int n = 0; n < nEn; n++) {
				ofile << setprecision(16) << SE[n] << '\t' << H[n] << endl;
			}
			ofile << endl << endl;
		} while (iter);
		////////////////
		ofile.close();

	}
	void RunMetroplisSE_MulRCE(double un0, double beta0, double alpha,
			double ENend, int iter) {
		double U0 = un0 * nN;
		Arg_io[AI_TEMPER] = beta0;
		for (int n = 0; n < nEn; n++) {
			SE[n] = E[n] * Arg_io[AI_TEMPER] + alpha / nN * ((E[n] - U0)
					* (E[n] - U0)) / 2.0;
		}
		SetSefun(SE, nEn);

		char filename[120] = { 0 };
		sprintf(fileprefix, "mulcep%d_%d_%d~%d~%lf_%lf_%lf_%lf_%d", nDim0,
				nDim1, int(cQ), nInitconf, un0, Arg_io[AI_TEMPER], alpha,
				ENend, iter);
		sprintf(filename, "%sStat.txt", fileprefix);
		std::ofstream ofile(filename);

		Gauge();
		unsigned long Samples = 10000;

		MetroplisTrialSe(Samples / 10);
		for (unsigned long s = 0; s < Samples; s++) {
			MetroplisTrialSe(nN);
			H[nEn - 1 + nCurHami] += 1;
		}
		for (int n = 0; n < nEn; n++) {
			ofile << setprecision(16) << SE[n] << '\t' << H[n] << endl;
		}
		ofile << endl << endl;

		double Eend = long(ENend * nN);
		double de = 0;
		double mBeta, b0;
		int EminN = -1, EmaxN, EmaxN2 = 0;

		for (int n = 0; n < nEn; n++) {
			de += H[n] * E[n];
			if (EminN < 0 && H[n] >= HIST_CUT) {
				EminN = n;
			}
		}
		de /= Samples;
		EmaxN = nEn - 1 + int(de);

		do {
			if (Eend >= E[EminN]) {
				if (H[EmaxN] >= HIST_CUT)
					iter--;
				EminN = nEn - 1 + int(Eend);
			}
			for (int n = EmaxN; n >= EminN; n--) {
				if (H[n] >= HIST_CUT) {
					EmaxN2 = n;
					break;
				}
			}
			for (int n = EminN; n <= EmaxN2; n++) {
				if (H[n] >= HIST_CUT) {
					SE[n] += log(H[n] / H[EmaxN2]);
				}
			}
			int k = 1;
			while (H[EminN + k] < HIST_CUT && k < EmaxN - EminN)
				k++;
			mBeta = (SE[EminN + k] - SE[EminN]) / k;

			for (int n = EminN; n > 0; n--) {
				b0 = mBeta + alpha / nN * (E[n] - E[EminN]);
				if (b0 < 0)
					b0 = 0;
				SE[n - 1] = SE[n] - b0;
			}

			for (int n = 0; n < nEn; n++) {
				H[n] = 0;
			}

			Samples = (EmaxN - EminN) * HIST_CUT * nDim1;
			cout << E[EminN] << '\t' << E[EmaxN2] << '\t' << E[EmaxN] << '\t'
					<< Samples << '\t' << iter << endl;

			MetroplisTrialSe(Samples / 10);
			for (unsigned long s = 0; s < Samples; s++) {
				MetroplisTrialSe(nN);
				H[nEn - 1 + nCurHami] += 1;
			}

			for (int n = 0; n < EmaxN; n++) {
				if (H[n] >= HIST_CUT) {
					EminN = n;
					break;
				}
			}
			SmoothH(EminN, EmaxN, 4);

			for (int n = 0; n < nEn; n++) {
				ofile << setprecision(16) << SE[n] << '\t' << H[n] << endl;
			}
			ofile << endl << endl;
		} while (iter);
		////////////////
		ofile.close();

	}
	void RunMetroplisSE_MulCE(double beta0, double ENend, int iter) {
		Arg_io[AI_TEMPER] = beta0;
		for (int n = 0; n < nEn; n++) {
			SE[n] = E[n] / Arg_io[AI_TEMPER];
		}
		SetSefun(SE, nEn);

		char filename[120] = { 0 };
		sprintf(fileprefix, "mulcep%d_%d_%d~%d~%lf_%lf_%d", nDim0, nDim1,
				int(cQ), nInitconf, Arg_io[AI_TEMPER], ENend, iter);
		sprintf(filename, "%sStat.txt", fileprefix);
		std::ofstream ofile(filename);

		Gauge();
		unsigned long Samples = 10000;

		MetroplisTrialSe(Samples / 10);
		for (unsigned long s = 0; s < Samples; s++) {
			MetroplisTrialSe(nN);
			H[nEn - 1 + nCurHami] += 1;
		}
		for (int n = 0; n < nEn; n++) {
			ofile << setprecision(16) << SE[n] << '\t' << H[n] << endl;
		}
		ofile << endl << endl;

		double Eend = long(ENend * nN);
		double de = 0;
		double dS;
		int EminN = -1, EmaxN, EmaxN2 = 0;

		for (int n = 0; n < nEn; n++) {
			de += H[n] * E[n];
			if (EminN < 0 && H[n] >= HIST_CUT) {
				EminN = n;
			}
		}
		de /= Samples;
		EmaxN = nEn - 1 + int(de);

		do {
			if (Eend >= E[EminN]) {
				if (H[EmaxN] >= HIST_CUT)
					iter--;
				EminN = nEn - 1 + int(Eend);
			}
			for (int n = EmaxN; n >= EminN; n--) {
				if (H[n] >= HIST_CUT) {
					EmaxN2 = n;
					break;
				}
			}
			for (int n = EminN; n <= EmaxN2; n++) {
				if (H[n] >= HIST_CUT) {
					SE[n] += log(H[n] / H[EmaxN2]);
				}
			}
			int k = 1;
			while (H[EminN + k] < HIST_CUT && k < EmaxN - EminN)
				k++;
			dS = (SE[EminN + k] - SE[EminN]) / k;

			for (int n = EminN; n > 0; n--) {
				SE[n - 1] = SE[n] - dS;
			}

			for (int n = 0; n < nEn; n++) {
				H[n] = 0;
			}

			Samples = (EmaxN - EminN) * HIST_CUT * nDim1;
			cout << E[EminN] << '\t' << E[EmaxN2] << '\t' << E[EmaxN] << '\t'
					<< Samples << '\t' << iter << endl;

			MetroplisTrialSe(Samples / 10);
			for (unsigned long s = 0; s < Samples; s++) {
				MetroplisTrialSe(nN);
				H[nEn - 1 + nCurHami] += 1;
			}

			for (int n = 0; n < EmaxN; n++) {
				if (H[n] >= HIST_CUT) {
					EminN = n;
					break;
				}
			}
			SmoothH(EminN, EmaxN, 4);
			for (int n = 0; n < nEn; n++) {
				ofile << setprecision(16) << SE[n] << '\t' << H[n] << endl;
			}
			ofile << endl << endl;
		} while (iter);
		////////////////
		ofile.close();

	}

	int RunStatQuan(char * file, unsigned long Samples, unsigned long rehead) {
		return 1;
	}

	~ CMulCEPotts() {
		delete[] States_p;
		delete[] SE;
		delete[] E;
		delete[] H;
	}
};

int main(int argc, char *argv[]) {
	clock_t __start = clock();

	if (argc == 8) {
		CMulCEPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
				atoi(argv[4]));
		p.RunMetroplisSE_MulCE(atof(argv[5]), atof(argv[6]), atoi(argv[7]));
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.fileprefix << std::endl;
	} else if (argc == 10) {
		CMulCEPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
				atoi(argv[4]));
		p.RunMetroplisSE_MulRCE(atof(argv[5]), atof(argv[6]), atof(argv[7]),
				atof(argv[8]), atoi(argv[9]));
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.fileprefix << std::endl;
	} else {
		cout << "<MulCEPotts> l0,l1,Q,initconf,B0,ENend,iter" << endl;//sampling at first time
		cout << "<MulCEPotts> l0,l1,Q,initconf,UN0,B0,Alpha,ENend,iter" << endl;//sampling at first time
		cout
				<< "<MulCEPotts> mulcep%d_%d_%d~%d~%lf_%lf_%dStatis.txt samples,stride"
				<< endl;//start from last time
		return -1;
	}
	return 0;
}

