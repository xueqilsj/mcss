/*
 ============================================================================
 Name        : NVTIsing.cpp
 Description : Sampling in NVT Ising model with Metroplis Algorithm
 ============================================================================
 */

#include "Ising.h"
#include "nvector.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;

class CNVTIsing: public Ising {
	double *E;

public:
	char filename[128];
	CNVTIsing(int dim0, int dim1, int initconf = 0) {
		InitParameters();
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(initconf);
		E = new double[nEn];
		for (int n = 0; n < nEn; n++) {
			E[n] = 4 * n - 2 * nN;
		}
	}
	~ CNVTIsing() {
		delete[] States_p;
		delete[] E;
	}
	void RunMetroplisTrial(double temp, unsigned long Samples,
			unsigned int Stride) {
		SetTemperature(temp);
		Gauge();
		sprintf(filename, "nvti%d_%d~%d~%lg_%ld_%dTraj.txt", nDim0, nDim1,
				nInitconf, temp, Samples, Stride);

		Stride *= nN;
		nvector<double> TR(Samples);

		for (unsigned long CurSample = 0; CurSample < Samples; CurSample++) {
			MetroplisTrial(Stride);

			TR.data[CurSample] = Arg_io[AI_HAMI];
		}
		TR.dump(filename);//save sample per 1e6

		nvector<char> cons(States_p, nN);
		sprintf(filename, "nvti%d_%d~%d~%lg_%ld_%dLastC.txt", nDim0, nDim1,
				nInitconf, temp, Samples, Stride / nN);//the last conf.
		cons.dump_bin(filename);
	}
	void RunHeatBathTrial(double temp, unsigned long Samples,
			unsigned int Stride) {
		SetHeatBath(temp);
		Gauge();

		sprintf(filename, "nvti%d_%d~%d~%lg_%ld_%dQuanFull.txt", nDim0, nDim1,
				nInitconf, temp, Samples, Stride);

		Stride *= nN;
		std::ofstream ofile(filename);
		ofile << setprecision(16);
		for (unsigned long CurSample = 0; CurSample < Samples; CurSample++) {
			HeatBathTrial(Stride);
			ofile << Arg_io[AI_HAMI] << endl;//<< '\t' << m2 << endl;
		}
		ofile.close();
		nvector<char> cons(States_p, nN);
		sprintf(filename, "nvti%d_%d~%d~%lg_%ld_%dLastC.txt", nDim0, nDim1,
				nInitconf, temp, Samples, Stride / nN);//the last conf.
		cons.dump_bin(filename);
	}

	void RunWolffCluster(double temp, unsigned long Samples,
			unsigned int WFsteps) {
		SetTemperature(temp);
		Gauge();
		ClusterStack_p = new int[nN];
		dPadd = 1 - exp(-2.0 / Arg_io[AI_TEMPER]);

		sprintf(filename, "nvti%d_%d~%d~%lg_%ld_%dTraj.txt", nDim0, nDim1,
				nInitconf, temp, Samples, WFsteps);

		nvector<double> TR(Samples);

		for (unsigned long CurSample = 0; CurSample < Samples; CurSample++) {
			WolffCluster(WFsteps);
			Gauge();
			TR.data[CurSample] = Arg_io[AI_HAMI];
		}
		TR.dump(filename);//save sample per 1e6

		nvector<char> cons(States_p, nN);
		sprintf(filename, "nvti%d_%d~%d~%lg_%ld_%dLastC.txt", nDim0, nDim1,
				nInitconf, temp, Samples, WFsteps);//the last conf.
		cons.dump_bin(filename);
		delete ClusterStack_p;
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

};
int main(int argc, char *argv[]) {
	if (argc == 8) {
		clock_t __start = clock();
		CNVTIsing p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
		if (argv[7][0] == 'm') {
			p.RunMetroplisTrial(atof(argv[4]), atoi(argv[5]), atoi(argv[6]));
		} else if (argv[7][0] == 'h') {
			p.RunHeatBathTrial(atof(argv[4]), atoi(argv[5]), atoi(argv[6]));
		} else if (argv[7][0] == 'w') {
			p.RunWolffCluster(atof(argv[4]), atoi(argv[5]), atoi(argv[6]));
		} else {
			cout << "unknow flag: " << argv[7] << endl;
		}
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;
	} else if (argc == 3) {
		int L0, L1, initconf;//!!!
		unsigned long Samples;
		unsigned int stride;
		double temp;
		int sv = sscanf(argv[1], "nvti%d_%d~%d~%lg_%ld_%dTraj", &L0, &L1,
				&initconf, &temp, &Samples, &stride);
		if (sv != 6) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CNVTIsing p(L0, L1, initconf);
		p.Arg_io[p.AI_TEMPER] = temp;
		p.RunStatQuan(argv[1], Samples, atol(argv[2]));

	} else {
		cout
				<< "<NVTIsing> dim0 dim0 initConf temperature samples stride m/h/w"
				<< endl;//stride =MCSweep/site
		cout << "<NVTIsing> nvti%d_%d~%d~%lg_%ld_%dTraj*.txt nLags" << endl;//nLags for AutoCorr
		return -1;
	}
	return 0;
}

