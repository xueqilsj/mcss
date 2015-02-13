#include <iostream>
#include "nvector.h"
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;
int main(int argc, char *argv[]) {
	if (argc == 6) {
		clock_t __start = clock();
		nvector<double> rdb;
		rdb.load(argv[1], atoi(argv[2]));
		unsigned long start = atol(argv[3]);
		unsigned long nblocks = atol(argv[4]);
		long tauMax = atol(argv[5]);
		//double mean, std, meanErr;
		//rdb.mean_std_err(&mean, &std, &meanErr, start);
		//cout << "%mean_std_err:" << mean << '\t' << std << '\t' << meanErr<< endl;

		double *tauCorr;
		double *blkAve;

		if (start < 0) {
			start = rdb.len / 30;
			if (start < 2) {
				cout << "Samples too little" << rdb.len << endl;
				return -1;
			}
			start = rdb.get_eq(start, 0);
		}

		long start2 = rdb.jack_knife_mean_error(nblocks, tauMax, tauCorr,
				blkAve, start);
		cout << "%jack_knife: start_mean_error:" << start2 << '\t'
				<< blkAve[nblocks] << '\t' << blkAve[nblocks + 1] << endl;
		nvector<double> rdb_e(tauCorr, (nblocks + 2) * tauMax);
		double times, error_time;
		long tau = rdb_e.integrated_time(nblocks, tauMax, &times, &error_time);
		cout << "%jack_knife: rho=" << times << "\trho_err=" << error_time
				<< "\ttau=" << tau << endl;

		cout << "jkTauCor=[";
		for (int t = 0; t < tauMax; t++)
			cout << tauCorr[nblocks * tauMax + t] << ' ';
		cout << "];" << endl;
		cout << "jkTauCor_ERR=[";
		for (int t = 0; t < tauMax; t++)
			cout << tauCorr[(nblocks + 1) * tauMax + t] << ' ';
		cout << "];" << endl;

		nvector<double> timeA(tauMax);
		nvector<double> errorA(tauMax);
		nvector<double> timeB(tauMax);
		nvector<double> errorB(tauMax);

		rdb_e.expo_times(nblocks, tauMax, timeA.data, errorA.data, timeB.data,
				errorB.data);
		cout << "A=[" << timeA << "];" << endl;
		cout << "AE=[" << errorA << "];" << endl;
		cout << "B=[" << timeB << "];" << endl;
		cout << "BE=[" << errorB << "];" << endl;
		//tauCorr[(nblo + 2) * tauMax]: nblo*tauMax,mean_Tau,error_Tau
		//blkAve[nblo+1]: (nblo-1)*nblo_1_ave,nblo_ave
		cout
				<< "figure;errorbar(jkTauCor,jkTauCor_ERR,'rs');ylabel('\\rho');set(gca,'yscale','log');\nfigure;errorbar(A(1:end-1),AE(1:end-1),'go');ylabel('A');\nfigure;errorbar(B(1:end-1),BE(1:end-1),'b^');ylabel('B');"
				<< endl;
		delete[] tauCorr;
		delete[] blkAve;
		cout << "%Escaped=" << ((double) (clock() - __start)) / CLOCKS_PER_SEC
				<< std::endl;
	} else {
		cout << "<AutoCorr> file col +/-start nblocks tauMax " << endl;//sampling at first time
	}
	return 1;
}
