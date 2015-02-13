/*
 ============================================================================
 Name        : NVEPotts.cpp
 Description : Sampling in NVE Potts model with Metroplis Algorithm
 ============================================================================
 */

#include "Potts.h"
#include "nvector.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

#define SAMPLES_ONCE 2000000

using namespace std;

class CNVEPotts: public Potts {
	double *E;
public:
	char filename[128];
	CNVEPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(initconf);
		E = new double[nEn];
		for (int n = 0; n < nEn; n++) {
			E[n] = nEn - n - 1;
			E[n] *= -1;
		}

	}
	void RunMetroplisCreutz(int nConsEnergy, unsigned long Samples,
			unsigned int Stride) {
		nEDemon = 2 * nN + nConsEnergy;
		if (nEDemon < 0 || nEDemon > 2 * nN) {
			cerr << "Constant Energy beyond Energy space at " << nConsEnergy
					<< endl;
			return;
		}
		Gauge();

		sprintf(filename, "nvep%d_%d_%d~%d~%d_%ld_%dTraj.txt", nDim0, nDim1,
				int(cQ), nInitconf, nConsEnergy, Samples, Stride);
		std::ofstream ofile(filename);
		Stride *= nN;
		/*
		 unsigned long nS = (Samples < SAMPLES_ONCE) ? Samples : SAMPLES_ONCE;

		 nvector<double> TR(nS);
		 nvector<int> EDemon(nS);
		 nvector<int> QNTR(nS * (cQ + 1));
		 nvector<double> QM2(nS);
		 unsigned int QDis[cQ + 1];
		 double m;
		 */
		for (unsigned long CurSample = 0; CurSample < Samples; CurSample++) {
			MetroplisTrialCreutzNVE(Stride);
			ofile << Arg_io[AI_HAMI] << '\t' << nEDemon << endl;
		}
		ofile.close();
	}
	void RunMetroplisTrialNVE(double RefU, unsigned long Samples,
			unsigned int Stride) {
		if (!SetRefU(RefU)) {
			cerr << "Failure in SetRefU()" << endl;
			return;
		}
		Gauge();

		sprintf(filename, "nvep%d_%d_%d~%d~%lf_%ld_%dTraj.txt", nDim0, nDim1,
				int(cQ), nInitconf, dRefU / nN, Samples, Stride);

		Stride *= nN;
		nvector<double> TR(Samples);
		nvector<double> MTR(Samples);

		for (unsigned long CurSample = 0; CurSample < Samples; CurSample++) {
			MetroplisTrialNVE(Stride);
			/*cout << nAccepted << '\t' << Arg_io[AI_HAMI] << '\t'
			 << Arg_io[AI_DELT_HAMI] << '\t'
			 << Arg_io[AI_EXP_DELT_HAMI] << '\t'
			 << Arg_io[AI_RANDOM] << endl;
			 */
			TR.data[CurSample] = Arg_io[AI_HAMI];
			MTR.data[CurSample] = Arg_io[AI_MAG];
		}
		nv_dump(TR, MTR, filename);//save sample per 1e6
	}

	void RunMetroplisSENVE(double RefU, unsigned long Samples,
			unsigned int Stride) {
		if (!SetRefU(RefU)) {
			cerr << "Failure in SetRefU()" << endl;
			return;
		}
		Gauge();

		char fileprefix[80] = { 0 };
		sprintf(fileprefix, "nvep%d_%d_%d~%d~%.14lf_%ld_%d", nDim0, nDim1,
				int(cQ), nInitconf, RefU, Samples, Stride);//use input
		sprintf(filename, "%sTraj.txt", fileprefix);

		int uN = 2 * nN + int(floor(dRefU));
		nvector<double> SE(nEn);
		for (int n = 0; n <= uN; n++) {
			SE.data[n] = (2.0 - nN) / 2.0 * log(dRefU - E[n]);
			//cout << n << "=\t" << SE.data[n] << endl;
		}
		double beta1 = (nN - 2.0) / (2.0 * (dRefU - E[uN]));
		for (int n = uN + 1; n < nEn; n++) {
			SE.data[n] = SE.data[n - 1] + beta1;
			//cout << n << "v\t" << SE[n] << endl;
		}
		SetSefun(SE.data, nEn);
		//for (int n = 0; n < nEn; n++)			cout << setprecision(16) << SE.data[n] << endl;

		std::ofstream ofile(filename);
		Stride *= nN;
		unsigned int *QDis = new unsigned int[cQ + 1];
		double m1, m2;
		nvector<char> cons(States_p, nN);
		unsigned long duS = Samples / 30;
		for (unsigned long s = 0; s < Samples; s++) {
			MetroplisTrialSe(Stride);
			m2 = QDistri(QDis, cQ + 1);
			m1 = cQ * QDis[0];
			m1 = (m1 / nN - 1) / (cQ - 1);
			ofile << setprecision(16) << Arg_io[AI_HAMI] << '\t'
					<< setprecision(6) << m1 << '\t' << m2 << endl;
			if (s % duS == 0) {
				sprintf(filename, "%sMidC%03ld_%lf_%lg_%lg.txt", fileprefix, s
						/ duS, Arg_io[AI_HAMI], m1, m2);//the last conf.
				cons.dump_bin(filename);
			}
		}
		delete[] QDis;
		ofile.close();
	}
	int RunStatQuanRay(char * file, unsigned long Samples, unsigned long nLags) {
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
		double un0 = dRefU;
		dRefU *= nN;
		double beta = (nN - 2.0) / (2.0 * (dRefU - mean));
		double alpha = beta * nN / (dRefU - mean);//Tc=1.4260624
		cout
				<< "#nInit\tUN0\tBeta\tAlpha\tWin\tN\tTeq\tnLags\tTau\tNeff\tmeanEd\tstdEd\tmeanE\tstdE"
				<< endl << nInitconf << '\t' << setprecision(16) << un0 << '\t'
				<< beta << '\t' << alpha << '\t' << win << '\t' << Esamples.len
				<< '\t' << uTeq << '\t' << uLags << '\t' << uTau << '\t'
				<< neff << '\t' << mean << '\t' << std << '\t' << mean1 << '\t'
				<< std1 << endl;
		return 0;
	}
	int RunStatQuan(char * file, unsigned long Samples, unsigned long nLags) {
		nvector<double> Esamples(Samples);
		nvector<double> EDemon(Samples);
		std::ifstream ifile(file);
		char eatline[256];
		for (unsigned int i = 0; i < EDemon.len; i++) {
			ifile >> Esamples.data[i] >> EDemon.data[i];
			ifile.getline(eatline, 256);
			if (ifile.fail()) {
				cerr << "#Failure in reading samples at" << i << endl;
				return -1;
			}
		}
		ifile.close();

		unsigned int uTau;//Time of AutoCorrelation (MCS/site)
		unsigned int uTeq;//Equilibrium at Teq
		unsigned int win = (unsigned int) (8 * pow(nN, 0.52));
		unsigned long uLags = EDemon.statis(nLags, win, uTeq, uTau);//~3*tau);
		if (uLags == 0)
			uTau = 1;
		double mean, std, mean1, std1;
		unsigned int neff = EDemon.mean_std(&mean, &std, false, uTeq, uTau);
		Esamples.mean_std(&mean1, &std1, false, uTeq, 1);
		Arg_io[AI_TEMPER] = 1.0 / log(1 + 1.0 / mean);

		cout
				<< "#nInit\tTemp\tWin\tN\tTeq\tnLags\tTau\tNeff\tmeanEd\tstdEd\tmeanE\tstdE"
				<< endl << nInitconf << '\t' << setprecision(16)
				<< Arg_io[AI_TEMPER] << '\t' << win << '\t' << EDemon.len
				<< '\t' << uTeq << '\t' << uLags << '\t' << uTau << '\t'
				<< neff << '\t' << mean << '\t' << std << '\t' << mean1 << '\t'
				<< std1 << endl;
	}
	~ CNVEPotts() {
		delete[] States_p;
		delete[] E;
	}
};
int main(int argc, char *argv[]) {
	if (argc == 8) {
		clock_t __start = clock();
		CNVEPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
		p.RunMetroplisSENVE(atof(argv[5]), atoi(argv[6]), atoi(argv[7]));//usually initConf=-1
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;
	} else if (argc == 3) {
		int L0, L1, Q, initconf;//!!!
		unsigned long Samples;
		unsigned int stride;
		double refU;
		int sv = sscanf(argv[1], "nvep%d_%d_%d~%d~%lf_%ld_%dTraj", &L0, &L1,
				&Q, &initconf, &refU, &Samples, &stride);
		if (sv != 7) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CNVEPotts p(L0, L1, Q, initconf);
		//p.nEDemon = 2 * p.nN + refU;
		p.dRefU = refU;
		p.RunStatQuanRay(argv[1], Samples, atol(argv[2]));

	} else {
		cout << "<NVEPotts> dim0 dim0 Q initConf refU samples stride" << endl;//stride =MCSweep/site
		cout << "<NVEPotts> nvep%d_%d_%d~%d~%lg_%ld_%dTraj*.txt nLags" << endl;//nLags for AutoCorr
		return -1;
	}
	return 0;
}

