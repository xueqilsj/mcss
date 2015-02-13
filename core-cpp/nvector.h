/*
 ============================================================================
 Name        : numvector.h
 Description :
 ============================================================================
 */

#ifndef _NVECTOR_H
#define _NVECTOR_H

#include <string>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

template<class T>
bool nv_split(const std::string &input, T &res, int c,
		const std::string& delimiters) {
	std::string token = "";
	std::string::size_type lastPos = input.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	std::string::size_type pos = input.find_first_of(delimiters, lastPos);

	while (std::string::npos != pos || std::string::npos != lastPos) {
		// Found a token, add it to the vector.
		if (c-- == 0) {
			token = input.substr(lastPos, pos - lastPos);
			break;
		}
		// Skip delimiters.  Note the "not_of"
		lastPos = input.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = input.find_first_of(delimiters, lastPos);
	}
	if (token == "")
		return false;
	std::istringstream ss(token);
	ss >> res;
	return true;
}

template<class DataType>
class nvector {
	bool own;

private:
	nvector(const nvector & other);
	nvector & operator =(const nvector & other);
public:
	DataType *data;
	std::string delim;
	unsigned int len;
public:
	nvector() {
		len = 0;
		data = NULL;
		delim = ",";
		own = false;
	}
	nvector(unsigned int _len, bool _own = true) {
		len = _len;
		data = new DataType[len];
		delim = ",";
		own = _own;
	}
	nvector(DataType *_data, unsigned int _len) {
		len = _len;
		data = _data;
		delim = ",";
		own = false;
	}
	///////////////
	~nvector() {
		if (own)
			delete[] data;
	}
	void fill(DataType val) {
		for (int i = 0; i < len; i++)
			data[i] = val;
	}
	void dump(const std::string &file, char linesep = '\n', unsigned int end =
			0) {
		std::ofstream ofile(file.c_str());
		for (unsigned int i = 0; i < len - end; i++) {
			ofile << data[i];
			if (i < len - 1 - end)
				ofile << linesep;
		}
		ofile.close();
	}
	void dump_bin(const char *file, std::ios::openmode mode = std::ios::out
			| std::ios::binary) {//| std::ios::app
		std::ofstream ofile(file, mode);
		ofile.write(data, len);
		ofile.close();
	}
	int load(const std::string& filename, int col = 0, bool o = true,
			char linesep = '\n', const std::string &sep = "\t") {
		std::ifstream In_f(filename.c_str());
		if (In_f.fail()) {
			std::cerr << "Failure in opening file:" << filename << std::endl;
			return -1;
		}

		char szLine[1024];
		DataType res;
		if (o) {
			std::vector<DataType> vData;
			while (In_f.getline(szLine, 1024, linesep)) {
				if (nv_split(std::string(szLine), res, col, sep)) {
					vData.push_back(res);
				}
			}
			if (own)
				delete[] data;
			data = new DataType[vData.size()];
			own = true;
			for (len = 0; len < vData.size(); len++)
				data[len] = vData[len];

		} else {
			unsigned int l = 0;
			while (In_f.getline(szLine, 1024, linesep)) {
				if (nv_split(std::string(szLine), res, col, sep)) {
					if (l >= len) {
						std::cerr << "too much items:" << l << std::endl;
						break;
					}
					data[l] = res;
					l++;
				}
			}
			len = l;

		}
		In_f.close();
		return len;
	}
	int load_bin(const char *file, unsigned int size, long offset = 0,
			bool _own = true) {
		std::ifstream infile(file, std::ifstream::binary);
		if (infile.fail()) {
			std::cerr << "Failure in opening file:" << file << std::endl;
			return -1;
		}
		infile.seekg(offset);
		if (_own) {
			if (own && len != 0) {
				delete[] data;
			}
			data = new DataType[size];
			own = true;
			infile.read(data, size);
		} else {
			if (size != len) {
				std::cerr << "too much items:" << len << std::endl;
				return -1;
			} else {
				infile.read(data, size);
			}
		}
		if (infile.eof()) {
			len = strlen(data);
			std::cout << "size" << len << std::endl;
		} else {
			len = size;
		}
		infile.close();
		return len;
	}
	friend std::ostream &operator<<(std::ostream &os, nvector &nv) {
		for (unsigned int i = 0; i < nv.len; i++) {
			os << nv.data[i];
			if (i < nv.len - 1)
				os << nv.delim;
		}
		return os;
	}
	DataType sum() {
		DataType m = 0;
		for (unsigned int i = 0; i < len; i++)
			m += data[i];
		return m;
	}
	double mean() {
		double m = 0;
		for (unsigned int i = 0; i < len; i++)
			m += data[i];
		return m / len;
	}
	double std() {
		double m = 0;
		double m2 = 0;
		for (unsigned int i = 0; i < len; i++) {
			m += data[i];
			m2 += data[i] * data[i];
		}
		return sqrt(m2 / len - (m / len) * (m / len));
	}

	unsigned int mean_std(double *mean, double *std, bool flag,//false=std,true=1/N
			unsigned int start = 0, unsigned int step = 1) {
		DataType *data_p = data + start;
		unsigned int N = len - start;
		*mean = 0;
		double m2 = 0;
		for (unsigned int k = 0; k < N; k += step) {
			*mean += data_p[k];
			m2 += data_p[k] * data_p[k];
		}
		unsigned int n = (N - 1) / step + 1;
		*mean /= n;
		m2 = m2 / n - (*mean) * (*mean);
		if (!flag)
			m2 = m2 * n / (n - 1);//
		*std = sqrt(m2);
		return n;
	}
	unsigned int mean_std_err(double *mean, double *std, double *meanErr,
			unsigned int start = 0, unsigned int tau = 0) {
		DataType *data_p = data + start;
		unsigned int N = len - start;
		*mean = 0;
		double m2 = 0;
		for (unsigned int k = 0; k < N; k++) {
			*mean += data_p[k];
			m2 += data_p[k] * data_p[k];
		}
		*mean /= N;
		m2 = m2 / N - (*mean) * (*mean);
		*meanErr = sqrt(m2 * (1 + 2 * tau) / (N - 1));
		*std = sqrt(m2 * N / (N - 1));
		return N;
	}
	unsigned int mean_std_err_blk(double *mean, double *std, double *meanErr,
			unsigned int start = 0, unsigned int blk = 1) {
		DataType *data_p = data + start;
		unsigned int N = len - start;
		*mean = 0;
		double lmean, m2 = 0;
		unsigned int nBlk = N / blk;
		for (unsigned int b = 0; b < N; b += nBlk) {
			lmean = 0;
			for (unsigned int i = b * blk; i < (b + 1) * blk; i++) {
				lmean += data_p[i];
			}
			lmean /= blk;
			*mean += lmean;
			m2 += lmean * lmean;
		}
		*mean /= nBlk;
		m2 = m2 / nBlk - (*mean) * (*mean);
		*meanErr = sqrt(m2 / (nBlk - 1));
		*std = sqrt(m2 * nBlk / (nBlk - 1));
		return nBlk;
	}
	unsigned int get_eq(unsigned int EQwin, unsigned int start = 0) {
		//get T_eq,EQwin>Tau correlation time
		DataType *data_p = data + start;
		unsigned int N = len - start;

		if (EQwin > N)
			return 0;

		double E0, S0, ErrorE0;
		//mean_std(&E0, &S0, false, start + N - EQwin);
		//ErrorE0 = S0 / sqrt(EQwin);
		mean_std_err(&E0, &S0, &ErrorE0, start + N - EQwin);

		double m = 0;
		unsigned int s = 0;
		while (s < N - EQwin) {
			m = 0;
			for (unsigned int k = s; k < s + EQwin; k++) {
				m += data_p[k];
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
	bool auto_corr(double *Xt, unsigned int nLags, unsigned int &tau,
			unsigned int start) {
		DataType *data_p = data + start;
		unsigned int N = len - start;

		if (nLags > N - 1)
			return false;//lags 'nLags' must not exceed 'Series' length - 1.

		double avm = 0;
		for (unsigned int i = 0; i < N; i++)
			avm += data_p[i];
		avm /= N;

		tau = 0;
		double mit, mitk;
		for (unsigned int k = 0; k <= nLags; k++) {
			//Xt[k] = _summ2(m, N, k, avm) / N;
			//def	_summ2(data, N, k, avm):
			Xt[k] = 0;
			mit = 0;
			mitk = 0;

			for (unsigned int it = 0; it < N - k; it++) {
				Xt[k] += data_p[it] * data_p[it + k];
				mit += data_p[it];
				mitk += data_p[it + k];
			}

			Xt[k] /= N - k;
			Xt[k] -= mit / (N - k) * mitk / (N - k);

			if (tau == 0 && (Xt[k] / Xt[0] < 0.135335283236613))//0.36787944117144233)//)//1/exp(1),)//1 / exp(1):
			{
				tau = (k + 1) / 2;//!=0
				break;
			}
		}
		return true;//although tau==0;
	}
	unsigned long statis(unsigned long uInitLags, unsigned int Win,
			unsigned int &uTeq, unsigned int &utau) //start,stop,step:MCS/site [0,M-1]
	{
		if (Win > len)
			return 0;
		uTeq = get_eq(Win, 0);//Win>Tau
		utau = 0;
		unsigned long uLags = uInitLags;
		double * dXT_p;
		do {
			dXT_p = new double[uLags + 1];
			if (!auto_corr(dXT_p, uLags, utau, uTeq)) {
				std::cerr << "#Failure in auto_corr(),'nLags':" << uLags
						<< " must not exceed 'Series' length - 1:" << len
						- uTeq << std::endl;
				return 0;
			} else if (utau == 0) {
				uLags *= 2;//(unsigned long) (pow(nN, 0.25));
			}
			delete[] dXT_p;
		} while (utau == 0);
		return uLags;
	}
	unsigned long jk_statis(int start, unsigned long nblo, long tauMax,
			unsigned int &uTeq, double &time, double &error_time,
			unsigned int &utau) //start,stop,step:MCS/site [0,M-1]
	{
		/*
		 if (start > len or en - start)
		 return 0;
		 uTeq = get_eq(Win, 0);//Win>Tau
		 utau = 0;
		 unsigned long uLags = uInitLags;
		 do {
		 //delete[] dXT_p;dXT_p = new double[uLags + 1];
		 double dXT_p[uLags + 1];
		 if (!auto_corr(dXT_p, uLags, utau, uTeq)) {
		 std::cerr << "#Failure in auto_corr(),'nLags':" << uLags
		 << " must not exceed 'Series' length - 1:" << len
		 - uTeq << std::endl;
		 return 0;
		 } else if (utau == 0) {
		 uLags *= 2;//(unsigned long) (pow(nN, 0.25));
		 }
		 } while (utau == 0);
		 return uLags;

		 double *tauCorr;
		 double *blkAve;

		 unsigned int win;
		 if (start > 0) {
		 win = (len - start) / 30;
		 get_eq(win, start);
		 } else {
		 win = len / 30;
		 get_eq(win, 0);
		 }

		 if (win < 2) {
		 cout << "Samples too little" << len << endl;
		 return -1;
		 }

		 start = rdb.get_eq(start, 0);

		 if (win_start < 0) {
		 start = len / 30;
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
		 */
		return 0;
	}

	//Calculate errors for a correlation function
	long jack_knife_mean_error(unsigned long nblo, long tauMax,
			double *&tauCorr, double*&blkAve, unsigned long start = 0) {
		//tauCorr[(nblo + 2) * tauMax]: nblo*tauMax,mean_Tau,Tau_err
		//blkAve[nblo+2]: (nblo-1)*nblo_1_ave,nblo_ave,nblo_ave_error

		unsigned long N = len - start;
		long block = N / nblo;
		N = block * nblo;
		DataType *data_p = data + (len - N);

		if (tauMax > block / 3) {
			std::cerr << "Failure tauMax too large:" << tauMax << ',' << block
					/ 3 << std::endl;
			return -1;
		}
		tauCorr = new double[(nblo + 2) * tauMax];
		blkAve = new double[nblo + 2];//average
		if (tauCorr == NULL) {
			std::cerr << "Failure: allocating memory for tauCorr:" << std::endl;
			return -1;
		}

		long tau, ti;
		unsigned long eb, e, b, bi;
		double sum, sum2;
		//average (N-1)  blocks
		for (b = 0; b < nblo; b++) {
			blkAve[b] = 0;
			for (e = 0; e < b * block; e++) {
				blkAve[b] += data_p[e];
			}
			for (e = (b + 1) * block; e < N; e++) {
				blkAve[b] += data_p[e];
			}

			blkAve[b] /= (nblo - 1) * block;

		}
		//average N  blocks
		/*
		 blkAve[nblo] = 0;
		 for (e = 0; e < N; e++) {
		 blkAve[nblo] += data_p[e];
		 }
		 blkAve[nblo] /= N;
		 */
		blkAve[nblo] = 0;
		for (b = 0; b < nblo; b++) {
			blkAve[nblo] += blkAve[b];
		}
		blkAve[nblo] /= nblo;

		//average error
		sum = 0;
		sum2 = 0;
		for (b = 0; b < nblo; b++) {
			sum += blkAve[b];
			sum2 += blkAve[b] * blkAve[b];
		}
		sum /= nblo;
		sum2 /= nblo;
		blkAve[nblo + 1] = sqrt((nblo - 1) * (sum2 - sum * sum));//save error to the last one

		//Corr
		for (b = 0; b < nblo; b++) {
			for (tau = 0; tau < tauMax; tau++) {
				eb = b * tauMax + tau;
				tauCorr[eb] = 0;
				//for every block every tau
				for (ti = 0; ti < block - tau; ti++) {
					for (bi = 0; bi < b; bi++) {
						tauCorr[eb] += (data_p[bi * block + ti] - blkAve[b])
								* (data_p[bi * block + ti + tau] - blkAve[b]);
					}
					for (bi = (b + 1); bi < nblo; bi++) {
						tauCorr[eb] += (data_p[bi * block + ti] - blkAve[b])
								* (data_p[bi * block + ti + tau] - blkAve[b]);
					}
				}
				tauCorr[eb] /= block - tau;
			}
		}

		//all Corr
		for (tau = 0; tau < tauMax; tau++) {
			eb = nblo * tauMax + tau;
			tauCorr[eb] = 0;
			//for every block every tau
			for (ti = 0; ti < block - tau; ti++) {
				for (bi = 0; bi < nblo; bi++) {
					tauCorr[eb] += (data_p[bi * block + ti] - blkAve[nblo])
							* (data_p[bi * block + ti + tau] - blkAve[nblo]);
				}
			}
			tauCorr[eb] /= block - tau;
			//std::cout << tauCorr[eb] << '\t';
		}

		for (tau = tauMax - 1; tau >= 0; tau--) //Go backwards to avoid spoiling t=0 value
		{//Go from C(t) to rho(t). Normalization is automatic.
			for (b = 0; b <= nblo; b++) {
				tauCorr[b * tauMax + tau] /= tauCorr[b * tauMax];
			}
		}

		//for error
		double pippo;
		for (tau = 0; tau < tauMax; tau++) {
			sum = sum2 = 0;
			for (b = 0; b < nblo; b++) {
				pippo = tauCorr[b * tauMax + tau];
				sum += pippo;
				sum2 += pippo * pippo;
			}
			sum /= nblo;
			sum2 /= nblo;
			tauCorr[(nblo + 1) * tauMax + tau] = sqrt((nblo - 1) * (sum2 - sum
					* sum));//save error to the last one
		}

		return len - N;//start

	}
	//Given an autocorrelation function and its error, this function
	//calculates its autocorrelation time with a window of size tauWin*tau_int
	//It return the time in variable "time" and the error estimate in "error_time"
	long integrated_time(unsigned long nblo, long tauMax, double *time,
			double *error_time, int tauWin = 6) {

		double sum, sum2, pippo, value;
		long tau;
		unsigned long b;
		sum = sum2 = 0;
		for (b = 0; b < nblo; b++) {
			pippo = 0.5;
			tau = 1;
			while ((tau < tauWin * pippo) && (tau < tauMax)) {
				value = data[b * tauMax + tau];
				if (fabs(value) > 2. * data[(nblo + 1) * tauMax + tau]) {
					pippo += value; //Add only if signal-to-noise ratio > 2
				}
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
		while ((tau < tauWin * pippo) && (tau < tauMax)) {
			value = data[nblo * tauMax + tau];
			if (fabs(value) > 2. * data[(nblo + 1) * tauMax + tau]) {
				pippo += value;
			}
			tau++;
		}
		if ((tau < tauWin * pippo) && (tau >= tauMax)) {
			std::cerr << "#warning: tau >= tauMax" << std::endl;
		}

		time[0] = pippo;
		return tau;
	}

	//Given an autocorrelation function and its error, this function
	//calculates  two estimators (and their errors) for the exponential
	//autocorrelation time
	void expo_times(unsigned long nblo, long tauMax, double *timeA,
			double *errorA, double *timeB, double *errorB) {
		//tauMax:timeA, errorA, timeB, errorB
		double sumA, sumA2, sumB, sumB2, pippo, value1, value2;
		long tau;
		unsigned long b;

		timeA[0] = 0;
		errorA[0] = 0;
		for (tau = 1; tau < tauMax - 1; tau++) {
			sumA = sumA2 = 0;
			for (b = 0; b < nblo; b++) {
				value1 = fabs(data[b * tauMax + tau]);//Take absolute values, in case
				pippo = -tau / log(value1);
				sumA += pippo;
				sumA2 += pippo * pippo;
			}
			value1 = fabs(data[nblo * tauMax + tau]);
			timeA[tau] = -tau / log(value1);
			sumA /= nblo;
			sumA2 /= nblo;
			errorA[tau] = sqrt((nblo - 1) * (sumA2 - sumA * sumA));
		}

		for (tau = 0; tau < tauMax - 1; tau++) {
			sumB = sumB2 = 0;
			for (b = 0; b < nblo; b++) {
				value1 = fabs(data[b * tauMax + tau]);//Take absolute values, in case
				value2 = fabs(data[b * tauMax + tau + 1]);//of anticorrelations.
				pippo = 1.0 / log(value1 / value2);
				sumB += pippo;
				sumB2 += pippo * pippo;
			}
			value1 = fabs(data[nblo * tauMax + tau]);
			value2 = fabs(data[nblo * tauMax + tau + 1]);
			timeB[tau] = 1.0 / log(value1 / value2);
			sumB /= nblo;
			sumB2 /= nblo;
			errorB[tau] = sqrt((nblo - 1) * (sumB2 - sumB * sumB));
			//std::cout << value1 << '\t' << value2 << '\t' << timeB[tau] << '\t'<< errorB[tau] << std::endl;
		}
	}

};

template<class T1, class T2>
int nv_dump(const nvector<T1> &v1, const nvector<T2> &v2, std::string file,
		char linesep = '\n', const std::string &sep = "\t") {
	if (v1.len != v2.len)
		return -1;
	std::ofstream ofile(file.c_str());
	for (unsigned int i = 0; i < v1.len; i++)
		ofile << v1.data[i] << sep << v2.data[i] << linesep;
	ofile.close();
	return v1.len;
}

template<class DataType>
void nv_dump(DataType * d, unsigned int l, std::string file, char linesep =
		'\n') {
	nvector<DataType> out(d, l);
	out.dump(file, linesep);
}
template<class DataType>
void nv_load(DataType * d, unsigned int l, std::string file, int col = 0,
		bool o = true, char linesep = '\n', const std::string &sep = "\t") {
	nvector<DataType> out(d, l);
	out.load(file, col, o, linesep, sep);
}

#endif
