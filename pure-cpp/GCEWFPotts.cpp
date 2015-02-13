/*
 ============================================================================
 Name        : GWFPotts.cpp
 Description : Sampling in Potts model with Wollf Cluster method
 ============================================================================
 */

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <time.h>

#include "Potts.h"
#include "nvector.h"

using namespace std;

class CGWFPotts: public Potts {
public:
	char filename[128];
	CGWFPotts() {
	}

	CGWFPotts(int dim0, int dim1, char Q, int initconf = 0) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(initconf);

	}
	void RunWolffClusterGCE(unsigned long Samples, unsigned long WFsteps,
			long nConfs) {
		sprintf(filename,
				"gcewfp%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ldTraj.txt",
				nDim0, nDim1, int(cQ), nInitconf, dBeta, dAlpha, dUn, Samples,
				WFsteps);
		std::ofstream ofile(filename);

		ClusterSeq_p = new int[nN];
		ClusterStack_p = new int[nN];

		Gauge();
		dGCEWF_AvgE = nCurHami;
		dGCEWF_SegAccumE = 0;
		nGCEWF_AvgC = 2000;
		//ran2.Srand(12387);
		//int(200 * pow(nN, 0.3));//      1024=12800, 16=1055
#ifdef _DEBUG
		cout << "nGCEWF_AvgC="<<nGCEWF_AvgC << endl;
#endif
		unsigned long snaps, uAvgC2 = 0;
		if (nConfs <= 0)
			snaps = Samples + 1;
		else
			snaps = Samples / nConfs;
		double Racc = 0;
		ofile << setprecision(16);
		for (unsigned long CurSample = 1; CurSample <= Samples; CurSample++) {
			WolffClusterGCE(WFsteps);
			if (CurSample % snaps == 0) {
				sprintf(
						filename,
						"gcewfp%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ld~%03ld_%d.conf",
						nDim0, nDim1, int(cQ), nInitconf, dBeta, dAlpha, dUn,
						Samples, WFsteps, CurSample / snaps, nCurHami);
				DumpConf(filename);

#ifndef _DEBUG
				dAccumCS /= (uCurMove / nGCEWF_AvgC - uAvgC2);
				uAvgC2 = uCurMove / nGCEWF_AvgC;
				Racc = exp(-dAlpha * nDeltaE * (nCurHami + nDeltaE / 2.0
						- dGCEWF_AvgE) / nN);
				cout << uCurMove << '\t' << csp << '\t' << dPadd << '\t'
						<< Racc << '\t' << nDeltaE << '\t' << dGCEWF_AvgE
						<< '\t' << dAccumCS << '\t' << nN / (((Racc > 1) ? 1
						: Racc) * dAccumCS) << endl;//2.0!!
				dAccumCS = 0;
#endif

			}
			ofile << nCurHami << endl;

		}
		ofile.close();
		delete[] ClusterSeq_p;
		delete[] ClusterStack_p;

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
		unsigned int win = Samples / 50;//(unsigned int) (8 * pow(nN, 0.52));
		unsigned long uLags = Esamples.statis(500, win, uTeq, uTau);//~3*tau);
		uTdrop = rehead;//Esamples.len / 10;//remove head 10%
		if (uTdrop < uTeq)
			uTdrop = uTeq;

		if (uLags == 0)
			uTau = 1;
		double mean, std, mean1, std1, meanErr;
		Esamples.mean_std_err(&mean, &std, &meanErr, uTdrop, uTau);
		double Beta = dBeta + dAlpha * (mean / nN - dUn);
		unsigned int neff = Esamples.mean_std(&mean1, &std1, false, uTdrop,
				uTau);

		cout << "#nInit\t" << setw(9) << "UN0"
				<< "\tB0\tAlpha\tWin\tN\tTeq\tnLags\tTau\tNeff\tmean\tstd\tBeta\tMerr\tmean1\tstd1"
				<< endl << nInitconf << '\t' << setw(9) << setprecision(16)
				<< dUn << '\t' << dBeta << '\t' << dAlpha << '\t' << win
				<< '\t' << Esamples.len << '\t' << uTdrop << '\t' << uLags
				<< '\t' << uTau << '\t' << neff << '\t' << mean << '\t' << std
				<< '\t' << Beta << '\t' << meanErr << '\t' << mean1 << '\t'
				<< std1 << '\t' << endl;
		return 0;
	}
	~ CGWFPotts() {
		delete[] States_p;
	}
};
int main(int argc, char *argv[]) {
	clock_t __start = clock();
	if (argc == 11) {
		CGWFPotts p(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
		p.dBeta = atof(argv[5]);
		p.dAlpha = atof(argv[6]);
		p.dUn = atof(argv[7]);
		p.RunWolffClusterGCE(atol(argv[8]), atol(argv[9]), atol(argv[10]));
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;
	} else if (argc == 8) {
		CGWFPotts p;
		if (p.FromConf(argv[1])) {
			p.dBeta = atof(argv[2]);
			p.dAlpha = atof(argv[3]);
			p.dUn = atof(argv[4]);
			p.RunWolffClusterGCE(atol(argv[5]), atol(argv[6]), atol(argv[7]));
			cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
					<< p.filename << std::endl;
			//delete States_p by ~CGWFPotts()
		}
	} else if (argc == 3) {
		int L0, L1, Q, initconf;//!!!
		unsigned long Samples;
		unsigned int stride;
		double beta, alpha, un;
		int sv = sscanf(argv[1], "gcewfp%d_%d_%d~%d~%lg_%lg_%lg_%ld_%dTraj",
				&L0, &L1, &Q, &initconf, &beta, &alpha, &un, &Samples, &stride);
		if (sv != 9) {
			cout << "Error input:" << argv[1] << '\t' << sv << endl;
			return -1;
		}
		CGWFPotts p(L0, L1, Q, initconf);
		p.dBeta = beta;
		p.dAlpha = alpha;
		p.dUn = un;
		p.RunStatQuan(argv[1], Samples, atol(argv[2]));

	} else {
		cout
				<< "<GCEWFPotts> dim0 dim0 Q initConf beta alpha un samples WFsteps nConf"
				<< endl;//WFstep :per Cluster
		cout
				<< "<GCEWFPotts> *%d_%d_%d~*Conf.txt B0 alpha dUn samples stride nConf"
				<< endl;//statistics
		cout
				<< "<GCEWFPotts> gcewfp%d_%d_%d~%d~%lg_%lg_%lg_%ld_%dTraj.txt nStart nBlock nTauMax"
				<< endl;//nLags for AutoCorr
		return -1;
	}
	return 0;

}

