/*
 ============================================================================
 Name        : WLXPotts.cpp
 Description : Sampling in Potts model with Wang-Landau method
 ============================================================================
 */

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <time.h>

#include "nvector.h"
#include "Potts.h"

using namespace std;

class CWLXPotts: public Potts {

public:
	char filename[128];
	char filePrefix[120];
	CWLXPotts(int dim0, int dim1, char Q) {
		InitParameters(Q);
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
		sprintf(filePrefix, "wlp%d_%d_%d~%lg_%lg_%d", nDim0, nDim1, int(cQ),
				lnfPrecision, flatRate, Stride);
		cout << "#WLXPotts\t" << nDim0 << '\t' << nDim1 << '\t' << int(cQ)
				<< '\t' << lnfPrecision << '\t' << flatRate << '\t' << Stride
				<< endl;
		unsigned long long samples = 0;
		clock_t __start = clock();
		Stride *= nN;
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
	void ContinueRunWL(double lastlnf, double lnfPrecision, double flatRate,
			unsigned int Stride = 1) {

		nvector<double> LnGe(nEn);
		nvector<unsigned long> GeN(nEn);
		SetWL(LnGe.data, GeN.data, nEn, lnfPrecision, flatRate);
		sprintf(filePrefix, "wlp%d_%d_%d~%lg_%lg_%d", nDim0, nDim1, int(cQ),
				lnfPrecision, flatRate, Stride);
		sprintf(filename, "%sLastC.txt", filePrefix);
		if (!LoadConf(filename)) {
			cout << "#Loading Conf. error" << filename << endl;
		}
		Gauge();
		Lnf = lastlnf;
		sprintf(filename, "%sLnGe_GeN.txt", filePrefix);
		ifstream ifile(filename);
		for (unsigned int i = 0; i < LnGe.len; i++)
			ifile >> LnGe.data[i] >> GeN.data[i];
		ifile.close();
		uCurMove = 0;
		for (int n = 0; n < nEn; n++) {
			GeN_p[n] = 0;
		}

		cout << "#WLXPotts\t" << nDim0 << '\t' << nDim1 << '\t' << int(cQ)
				<< '\t' << lnfPrecision << '\t' << flatRate << '\t' << Stride
				<< endl;
		unsigned long long samples = 0;
		clock_t __start = clock();
		Stride *= nN;
		while (Lnf > LnfPrecision) {
			WangLandauTrial(Stride);
			if (WLCheckBelowAVG()) {
				Lnf /= 2.0;
				samples += uCurMove;
				cout << Lnf << '\t' << uCurMove << '\t' << samples << '\t'
						<< ((double) (clock() - __start)) / CLOCKS_PER_SEC
						<< endl;

				sprintf(filename, "%sLnGe_GeN.txt", filePrefix);
				//nv_dump(LnGe, GeN, filename);

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
	void TestWLSE(char *fileLnGe, unsigned long Samples, unsigned int Stride) {
		InitConf(0);
		sprintf(filename, "%sLastC.txt", filePrefix);
		nvector<char> cons(States_p, nN);
		cons.load_bin(filename, nN, 0, false);//set last conf
		if (cons.len != nN) {
			cerr << "#Failure in loading" << filename << endl;
			return;
		}
		Gauge();
		nvector<double> vLnGe;
		vLnGe.load(fileLnGe, 0);
		if (!SetSefun(vLnGe.data, vLnGe.len)) {
			return;
		}

		unsigned int MCS = nN * Stride;
		nvector<double> TR(Samples);
		nvector<double> MTR(Samples);
		sprintf(filename, "%s~%ld_%dQuanFull.txt", filePrefix, Samples, Stride);
		unsigned long CurSample = 0, otSamples = 0;
		do {
			if (Samples - CurSample > 4000000)//4e6~2.2hours
				otSamples += 4000000;
			else
				otSamples = Samples;

			for (; CurSample < otSamples; CurSample++) {
				MetroplisTrialSe(MCS);

				TR.data[CurSample] = Arg_io[AI_HAMI];
				MTR.data[CurSample] = Arg_io[AI_MAG];
			}
			nv_dump(TR, MTR, filename);//save sample per 1e6
		} while (CurSample < Samples);

		sprintf(filename, "%s~%ld_%dLastC.txt", filePrefix, Samples, Stride);//the last conf.
		cons.dump_bin(filename);

	}
	~ CWLXPotts() {
		delete[] States_p;
	}
};

int main(int argc, char *argv[]) {
	int dim0, dim1, q;
	double lnfPrecision, flatRate;
	unsigned int stride;

	if (argc == 7) {//dim0 dim0 Q lnfPrecision,flatRate stride
		CWLXPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
		p.RunWL(atof(argv[4]), atof(argv[5]), atoi(argv[6]));
	} else if (argc == 8) {
		CWLXPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
		p.ContinueRunWL(atof(argv[7]), atof(argv[4]), atof(argv[5]), atoi(
				argv[6]));

	} else if (argc == 4) {//LnGefilename samples stride
		if (sscanf(argv[1], "wlp%d_%d_%d~%lg_%lg_%dLnGe_GeN.txt", &dim0, &dim1,
				&q, &lnfPrecision, &flatRate, &stride) != 6) {
			cout << "Unknown input:" << argv[1] << endl;
			return -1;
		}
		CWLXPotts p(dim0, dim1, q);
		sprintf(p.filePrefix, "wlp%d_%d_%d~%lg_%lg_%d", dim0, dim1, q,
				lnfPrecision, flatRate, stride);
		p.TestWLSE(argv[1], atoi(argv[2]), atoi(argv[3]));

	} else {
		cout << "<WLXPotts> dim0 dim1 Q lnfPrecision flatRate stride" << endl;
		cout << "<WLXPotts> dim0 dim1 Q lnfPrecision flatRate stride LastLnf"
				<< endl;
		cout << "<WLXPotts> LnGefilename samples stride" << endl;
		return -1;
	}

	return 0;
}
