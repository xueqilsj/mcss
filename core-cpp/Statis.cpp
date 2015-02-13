#include "Statis.h"
#include <cmath>
#include <iostream>

unsigned int MeanStd(double *TR_p, unsigned int N, double *mean, double *std,
		bool flag, unsigned int step) {
	*mean = 0;
	double m2 = 0;
	for (unsigned int k = 0; k < N; k += step) {
		*mean += TR_p[k];
		m2 += TR_p[k] * TR_p[k];
	}
	unsigned int n = (N - 1) / step + 1;
	*mean /= n;
	m2 = m2 / n - (*mean) * (*mean);
	if (!flag)
		m2 = m2 * n / (n - 1);//
	*std = sqrt(m2);
	return n;
}
double MeanStdErr(double *TR_p, unsigned int N, double *mean, double *std,
		double *meanErr, unsigned int tau) {
	*mean = 0;
	double m2 = 0;
	for (unsigned int k = 0; k < N; k++) {
		*mean += TR_p[k];
		m2 += TR_p[k] * TR_p[k];
	}
	*mean /= N;
	m2 = m2 / N - (*mean) * (*mean);
	*meanErr = sqrt(m2 * (1 + 2 * tau) / (N - 1));
	*std = sqrt(m2 * N / (N - 1));
	return m2;
}

int Jackknife_MeanStdErrTau_I(int *data_p, unsigned long len,
		unsigned int nblo, double *tauCorr, int tauMax, double *mean,
		double *std, double *meanErr, double *ptau) {
	//return tau
	if (tauCorr == NULL) {
		std::cerr << "Failure: allocating memory for tauCorr:" << std::endl;
		return -1;
	}

	unsigned long N = len;
	long block = N / nblo;
	N = block * nblo;
	data_p = data_p + (len - N);

	if (tauMax > block) {
		std::cerr << "Failure tauMax too large:" << tauMax << ',' << block
				<< std::endl;
		return -1;
	}
	int tau, ti;
	unsigned long e;

	//average N  blocks
	*mean = 0;
	double m2 = 0;
	for (e = 0; e < N; e++) {
		*mean += data_p[e];
		m2 += data_p[e] * data_p[e];
	}
	*mean /= N;
	m2 = m2 / N - (*mean) * (*mean);
	*std = sqrt(m2 * N / (N - 1));

	//all Corr
	for (tau = 0; tau < tauMax; tau++) {
		tauCorr[tau] = 0;
		//for every block every tau
		for (ti = 0; ti < block - tau; ti++) {
			for (unsigned long bi = 0; bi < nblo; bi++) {
				tauCorr[tau] += (data_p[bi * block + ti] - *mean) * (data_p[bi
						* block + ti + tau] - *mean);
			}
		}
		tauCorr[tau] /= block - tau;

	}
	for (tau = tauMax - 1; tau >= 0; tau--) //Go backwards to avoid spoiling t=0 value
	{//Go from C(t) to rho(t). Normalization is automatic.
		tauCorr[tau] /= tauCorr[0];
	}

	double pippo, value;
	long tauWin = 4;
	pippo = 0.5;
	tau = 1;
	while ((tau < tauWin * pippo) && (tau < tauMax)) {
		value = tauCorr[tau];
		if (fabs(value) > 0.36787944117144233) {
			pippo += value;
		}
		tau++;
	}
	if ((tau < tauWin * pippo) && (tau >= tauMax)) {
		std::cerr << "#warning: tau >= tauMax" << std::endl;
	}
	*meanErr = sqrt(m2 * (2 * pippo) / (N - 1));
	*ptau = pippo;
	return tau;
}
int Jackknife_MeanStdErrTau_D(double *data_p, unsigned long len,
		unsigned int nblo, double *tauCorr, int tauMax, double *mean,
		double *std, double *meanErr, double *ptau) {
	//return tau
	if (tauCorr == NULL) {
		std::cerr << "Failure: allocating memory for tauCorr:" << std::endl;
		return -1;
	}

	unsigned long N = len;
	long block = N / nblo;
	N = block * nblo;
	data_p = data_p + (len - N);

	if (tauMax > block) {
		std::cerr << "Failure tauMax too large:" << tauMax << ',' << block
				<< std::endl;
		return -1;
	}
	int tau, ti;
	unsigned long e;

	//average N  blocks
	*mean = 0;
	double m2 = 0;
	for (e = 0; e < N; e++) {
		*mean += data_p[e];
		m2 += data_p[e] * data_p[e];
	}
	*mean /= N;
	m2 = m2 / N - (*mean) * (*mean);
	*std = sqrt(m2 * N / (N - 1));

	//all Corr
	for (tau = 0; tau < tauMax; tau++) {
		tauCorr[tau] = 0;
		//for every block every tau
		for (ti = 0; ti < block - tau; ti++) {
			for (unsigned long bi = 0; bi < nblo; bi++) {
				tauCorr[tau] += (data_p[bi * block + ti] - *mean) * (data_p[bi
						* block + ti + tau] - *mean);
			}
		}
		tauCorr[tau] /= block - tau;

	}
	for (tau = tauMax - 1; tau >= 0; tau--) //Go backwards to avoid spoiling t=0 value
	{//Go from C(t) to rho(t). Normalization is automatic.
		tauCorr[tau] /= tauCorr[0];
	}

	double pippo, value;
	long tauWin = 4;
	pippo = 0.5;
	tau = 1;
	while ((tau < tauWin * pippo) && (tau < tauMax)) {
		value = tauCorr[tau];
		if (fabs(value) > 0.36787944117144233) {
			pippo += value;
		}
		tau++;
	}
	if ((tau < tauWin * pippo) && (tau >= tauMax)) {
		std::cerr << "#warning: tau >= tauMax" << std::endl;
	}
	*meanErr = sqrt(m2 * (2 * pippo) / (N - 1));
	*ptau = pippo;
	return tau;
}
double Jackknife_Xj(double *TR_p, unsigned int N, unsigned int j,
		unsigned int step) {//step=1
// X_i^J=[\sigma_{j!=i} x_j]/(N-1)
// f_i^J=f(x_i^J)
// \bar f^J= [\sigma f_i^J]/N
// \sigma_{fjerr}^2=(N-1)(\bar{f^j)^2}-(\bar f^J)^2 )
//j: start with 0, you can test the estimator of X^2

	unsigned int n = (N - 1) / step + 1;
	if (j >= n) {
		std::cout << "Jackknife_Xj:index override" << N << '\t' << j
				<< std::endl;
		return 0;
	}

	double m = 0;
	for (unsigned int k = 0; k < j; k += step) {
		m += TR_p[k];
	}
	for (unsigned int k = j + 1; k < N; k += step) {
		m += TR_p[k];
	}
	return m / (n - 1);
}
int AutoCor(double * TR_p, int N, double *Xt, int nLags) {
	if (nLags > N - 1)
		return -1;//lags 'nLags' must not exceed 'Series' length - 1.

	double avm = 0;
	for (int i = 0; i < N; i++)
		avm += TR_p[i];
	avm /= N;

	int ta = 0;
	double mit, mitk;
	for (int k = 0; k <= nLags; k++) {
		//Xt[k] = _summ2(m, N, k, avm) / N;
		//def	_summ2(TR, N, k, avm):
		Xt[k] = 0;
		mit = 0;
		mitk = 0;

		for (int it = 0; it < N - k; it++) {
			Xt[k] += TR_p[it] * TR_p[it + k];
			mit += TR_p[it];
			mitk += TR_p[it + k];
		}

		Xt[k] /= N - k;
		Xt[k] -= mit / (N - k) * mitk / (N - k);

		if (ta == 0 && Xt[k] / Xt[0] < 0.135335283236613)//0.36787944117144233)//1/exp(1),)//1 / exp(1):
		{
			ta = k / 2;
		}
	}
	return ta;
}

//Given an autocorrelation function and its error, this function
//calculates its autocorrelation time with a window of size WINDOW*tau_int
//It return the time in variable "time" and the error estimate in "error_time"
/*
 void integrated_time(double *corr, double *errorcorr, int nblo, int maxtau,
 double *time, double *error_time) {
 double sum, sum2, pippo, value;
 int myblo, tau;

 sum = sum2 = 0;
 for (myblo = 0; myblo < nblo; myblo++) {
 pippo = 0.5;
 tau = 1;
 while ((tau < WINDOW * pippo) && (tau < maxtau)) {
 value = corr[myblo * MAXTAU + tau];
 if (fabs(value) > 2. * errorcorr[tau])
 pippo += value; //Add only if signal-to-noise ratio > 2
 tau++;
 }
 sum += pippo;
 sum2 += pippo * pippo;
 }
 sum /= nblo;
 sum2 /= nblo;
 error_time[0] = sqrt((nblo - 1) * (sum2 - sum * sum));
 pippo = 0.5;
 tau = 1;
 while ((tau < WINDOW * pippo) && (tau < maxtau)) {
 value = corr[nblo * MAXTAU + tau];
 if (fabs(value) > 2. * errorcorr[tau])
 pippo += value;
 tau++;
 }
 time[0] = pippo;
 }
 */
//Given an autocorrelation function and its error, this function
//calculates  two estimators (and their errors) for the exponential
//autocorrelation time
/*
 void expo_times(double *corr, int nblo, int maxtau, double *timeA,
 double *errorA, double *timeB, double *errorB) {

 double sumA, sumA2, sumB, sumB2, pippo, value1, value2;
 int myblo, tau;
 for (tau = 1; tau < maxtau - 1; tau++) {
 sumA = sumB = sumA2 = sumB2 = 0;
 for (myblo = 0; myblo < nblo; myblo++) {
 value1 = fabs(corr[myblo * MAXTAU + tau]);//Take absolute values, in case
 value2 = fabs(corr[myblo * MAXTAU + tau + 1]);//of anticorrelations.
 pippo = -tau / log(value1);
 sumA += pippo;
 sumA2 += pippo * pippo;
 pippo = 1.0 / log(value1 / value2);
 sumB += pippo;
 sumB2 += pippo * pippo;
 }
 value1 = fabs(corr[nblo * MAXTAU + tau]);
 value2 = fabs(corr[nblo * MAXTAU + tau + 1]);
 timeA[tau] = -tau / log(value1);
 timeB[tau] = 1.0 / log(value1 / value2);
 sumA /= nblo;
 sumA2 /= nblo;
 sumB /= nblo;
 sumB2 /= nblo;
 errorA[tau] = sqrt((nblo - 1) * (sumA2 - sumA * sumA));
 errorB[tau] = sqrt((nblo - 1) * (sumB2 - sumB * sumB));
 }
 }
 */
unsigned int MoveAvg(double * TR, unsigned int N, unsigned int EQwin) {
	//get T_eq,EQwin>Tau correlation time
	if (EQwin > N)
		return 0;

	double E0, S0, ErrorE0;
	MeanStd(TR + N - EQwin, EQwin, &E0, &S0);
	ErrorE0 = S0 / sqrt(double(EQwin));

	double m = 0;
	unsigned int s = 0;
	while (s < N - EQwin) {
		m = 0;
		for (unsigned int k = s; k < s + EQwin; k++) {
			m += TR[k];
		}
		m /= EQwin;
		if (fabs(m - E0) < ErrorE0) {
			break;

		} else {
			s++;
		}
	}
	return s;
}
/*
 HisWeigh2NVT(doulbe t,double *E,double *SE,double *UList,double *FList,double * SList ,double * CList)
 {
 double e2_l = 0,e_l = 0,z_l = 0;
 l = SE[Level / 2] - E[Level / 2] / t
 for n in range (0, Level):
 x = SE[n] - E[n] / t - l
 try:
 p = exp(x)
 except OverflowError:
 print "error: exp(x),x=", x
 print "l=", l
 raise
 z_l += p
 e_l += E[n] * p
 e2_l += E[n] * E[n] * p;
 U = e_l / z_l
 F = -t * (log(z_l) + l)

 C = (e2_l / z_l - U * U) / (t * t)
 S = (U - F) / t
 UList.append(U / _nN)
 FList.append(F / _nN)
 CList.append(C / _nN)
 SList.append(S / _nN)
 */
