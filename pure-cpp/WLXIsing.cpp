/*
 ============================================================================
 Name        : WLXIsing.cpp
 Description : Sampling in Ising model with Wang-Landau method
 ============================================================================
 */

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <time.h>
#include "math.h"
#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))

#include "nvector.h"
#include "Ising.h"

using namespace std;

class CWLXIsing: public Ising {

public:
	char filename[128];
	char filePrefix[120];
	CWLXIsing(int dim0, int dim1) {
		InitParameters();
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
	}
	void RunWL(double lnfPrecision, double flatRate, unsigned int Stride = 1,
			int d = 0) {
		InitConf(d);
		Gauge();
		nvector<double> LnGe(nEn);
		nvector<unsigned long> GeN(nEn);
		SetWL(LnGe.data, GeN.data, nEn, lnfPrecision, flatRate);
		sprintf(filePrefix, "wli%d_%d~%lg_%lg_%d", nDim0, nDim1, lnfPrecision,
				flatRate, Stride);
		cout << "#WLXIsing\t" << nDim0 << '\t' << nDim1 << '\t' << lnfPrecision
				<< '\t' << flatRate << '\t' << Stride << endl;
		unsigned long long samples = 0;
		clock_t __start = clock();
		Stride *= nN;
		nLowEBound = nHighEBound = 2 * nN;

		WarmupToBound();
		while (Lnf > LnfPrecision) {
			WangLandauTrial(Stride);
			if (WLCheckBelowAVG()) {
				Lnf /= 2.0;
				samples += uCurMove;
				cout << Lnf << '\t' << uCurMove << '\t' << samples << '\t'
						<< ((double) (clock() - __start)) / CLOCKS_PER_SEC
						<< endl;
				sprintf(filename, "%sLnGe_GeN.txt", filePrefix);
				ofstream ofile(filename);
				for (unsigned int i = 0; i < LnGe.len; i++)
					ofile << setprecision(16) << LnGe.data[i] << '\t'
							<< GeN.data[i] << endl;
				ofile.close();

				sprintf(filename, "%sLastC.txt", filePrefix);
				DumpConf(filename);

				if (Lnf > LnfPrecision) {
					WLResetGeN();
				}
			}
		}

	}
	void RunWLBound(double lb, double ub, double lnfPrecision, double flatRate,
			unsigned int Stride = 1, int d = 0) {
		double eE = round((lb + 2) * nN / 4);
		nLowEBound = int(4 * eE - 2 * nN);
		eE = round((ub + 2) * nN / 4);
		nHighEBound = int(4 * eE - 2 * nN);

		InitConf(d);
		Gauge();
		WarmupToBound();

		nvector<double> LnGe(nEn);
		nvector<unsigned long> GeN(nEn);

		SetWL(LnGe.data, GeN.data, nEn, lnfPrecision, flatRate);
		sprintf(filePrefix, "wli%d_%d~%lg_%lg_%d~%d_%d", nDim0, nDim1,
				lnfPrecision, flatRate, Stride, nLowEBound, nHighEBound);
		cout << "#WLXIsing\t" << nDim0 << '\t' << nDim1 << '\t' << lnfPrecision
				<< '\t' << flatRate << '\t' << Stride << '\t' << nLowEBound
				<< '\t' << nHighEBound << endl;
		unsigned long long samples = 0;
		clock_t __start = clock();
		Stride *= nN;

		while (Lnf > LnfPrecision) {
			WangLandauTrialBound(Stride);
			if (WLCheckBelowAVGBound()) {
				Lnf /= 2.0;
				samples += uCurMove;
				cout << Lnf << '\t' << uCurMove << '\t' << samples << '\t'
						<< ((double) (clock() - __start)) / CLOCKS_PER_SEC
						<< endl;
				sprintf(filename, "%sLnGe_GeN.txt", filePrefix);
				ofstream ofile(filename);
				for (unsigned int i = 0; i < LnGe.len; i++)
					ofile << setprecision(16) << LnGe.data[i] << '\t'
							<< GeN.data[i] << endl;
				ofile.close();

				sprintf(filename, "%sLastC.txt", filePrefix);
				DumpConf(filename);

				if (Lnf > LnfPrecision) {
					WLResetGeNBound();
				}
			}
		}

	}
	~ CWLXIsing() {
		delete[] States_p;
	}
};

int main(int argc, char *argv[]) {
	//int dim0, dim1;
	//double lnfPrecision, flatRate;
	//unsigned int stride;

	if (argc == 6) {//dim0 dim0 Q lnfPrecision,flatRate stride
		CWLXIsing p(atoi(argv[1]), atoi(argv[2]));
		//p.RunWL(atof(argv[3]), atof(argv[4]), atoi(argv[5]));
		p.RunWLBound(-1.9, 1.9, atof(argv[3]), atof(argv[4]), atoi(argv[5]));
	} else {
		cout << "<WLXIsing> dim0 dim1 lnfPrecision flatRate stride" << endl;
		return -1;
	}

	return 0;
}
