#ifndef _STATIS_H_
#define _STATIS_H_

unsigned int MeanStd(double *TR_p, unsigned int N, double *mean, double *std,
		bool flag = false, unsigned int step = 1);//false=std,true=1/N
double MeanStdErr(double *TR_p, unsigned int N, double *mean, double *std,
		double *meanErr, unsigned int tau = 0);
double Jackknife_Xj(double *TR_p, unsigned int N, unsigned int j,
		unsigned int step = 1);//step=1,j count from 0;
unsigned int MoveAvg(double * TR, unsigned int N, unsigned int EQwin);
int AutoCor(double * TR_p, int N, double *Xt, int nLags);//Xt[nLags+1]
//can I resort to the keyword "export" ?
int Jackknife_MeanStdErrTau_I(int *data_p, unsigned long len,
		unsigned int nblo, double *tauCorr, int tauMax, double *mean,
		double *std, double *meanErr, double *ptau);
int Jackknife_MeanStdErrTau_D(double *data_p, unsigned long len,
		unsigned int nblo, double *tauCorr, int tauMax, double *mean,
		double *std, double *meanErr, double *ptau);
#endif
