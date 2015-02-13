#include <omp.h>
#include "Potts.h"
#include "nvector.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;

class TemperOMP_GCEPotts: public Potts {
protected:
	ofstream _fQuan;
	unsigned long _uSamples, _uCurSample, _uStride, _uSampStep;
	int nRank;
	///////////
public:
	char filename[128];
	unsigned long uAccepted;
	unsigned long uChecked;
	bool outAcc;
	TemperOMP_GCEPotts() {
	}
	virtual ~TemperOMP_GCEPotts() {
		delete[] States_p;
	}
	void _OpenOutPut(std::ios::openmode mode) {
		sprintf(filename,
				"togcep%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ldTraj.txt",
				nDim0, nDim1, int(cQ), nInitconf, dBeta, dAlpha, dUn,
				_uSamples, _uStride);
		cout << filename << endl;
		_fQuan.open(filename, mode);
	}
	void _CloseOutPut() {
		_fQuan.close();
	}
	///////////
	int Initiate(int dim0, int dim1, char Q, int initconf, double dBeta0,
			double alpha, double un0, unsigned long Samples,
			unsigned long Stride, int doflag, int Rank) {

		if (Stride == 0) {
			cout << "please set _uStride nonzero" << endl;
			return -1;
		}
		if (doflag == 0) {
			sprintf(filename,
					"togcep%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ldLastC.txt",
					dim0, dim1, int(Q), initconf, dBeta0, alpha, un0, Samples,
					Stride);
			if (!FromConf(filename)) {
				cerr << "#LoadConf Error:" << filename << endl;
				return -1;
			}
		} else {
			InitParameters(Q);
			States_p = new char[dim0 * dim1 + 1];
			SetConf(States_p, dim0, dim1);
			InitConf(initconf);
		}
		_uSamples = Samples;
		_uStride = Stride;
		_uSampStep = 0;
		_uCurSample = 0;
		nRank = Rank;
		uAccepted = 0;
		uChecked = 0;
		outAcc = false;
		dUn = un0;
		dBeta = dBeta0;
		dAlpha = alpha;

		time_t now = time(NULL);
		now /= 31;
		now *= 18653323 + nRank * 1231532129; //< Let every replica different random seed
		ran.Srand(now);
		return 1;
	}
	int RunStatQuan(int dim0, int dim1, char Q, int initconf, double dBeta0,
			double alpha, double un0, unsigned long Samples,
			unsigned long Stride, unsigned long rehead, int nRank) {

		if (Stride == 0) {
			cout << "please set _uStride nonzero" << endl;
			return -1;
		}
		_uSamples = Samples;
		_uStride = Stride;
		nDim0 = dim0;
		nDim1 = dim1;
		cQ = Q;
		nInitconf = initconf;
		dUn = un0;
		dBeta = dBeta0;
		dAlpha = alpha;
		sprintf(filename,
				"togcep%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ldTraj.txt",
				nDim0, nDim1, int(cQ), nInitconf, dBeta, dAlpha, dUn,
				_uSamples, _uStride);

		nvector<double> Esamples;
		Esamples.load(filename, 0);
		if (Esamples.len != Samples) {
			cerr << "#Reading Error:" << filename << endl;
			return -1;
		}
		unsigned int uTau;//Time of AutoCorrelation (MCS/site)
		unsigned int uTeq, uTdrop;//Equilibrium at Teq
		unsigned int win = (unsigned int) (8 * pow(nN, 0.52));
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

		sprintf(filename,
				"togcep%d_%d_%d~%d~%.15lg_%.15lg_%02d_%ld_%ldStat.txt", nDim0,
				nDim1, int(cQ), nInitconf, dBeta, dAlpha, nRank, _uSamples,
				_uStride);
		ofstream _fStatOut(filename);
		_fStatOut << "#nInit\t" << setw(9) << "dUn"
				<< "\tB0\tdAlpha\tWin\tN\tTeq\tnLags\tTau\tNeff\tmean\tstd\tBeta\tMerr\tmean1\tstd1"
				<< endl << nInitconf << '\t' << setw(9) << setprecision(16)
				<< dUn << '\t' << dBeta << '\t' << dAlpha << '\t' << win
				<< '\t' << Esamples.len << '\t' << uTdrop << '\t' << uLags
				<< '\t' << uTau << '\t' << neff << '\t' << mean << '\t' << std
				<< '\t' << Beta << '\t' << meanErr << '\t' << mean1 << '\t'
				<< std1 << endl;
		_fStatOut.close();
		return 1;
	}
	virtual bool RunMetroplisTrialGCE(unsigned long _uOneTimeMoves) {
		bool bContinue = true;
		unsigned int QDis[cQ + 1];
		double m1, m2;
		_fQuan << setprecision(16);
		while ((_uSampStep + _uOneTimeMoves) >= nN * _uStride) {//sampling
			MetroplisTrialGCE(nN * _uStride - _uSampStep);
			_uOneTimeMoves -= (nN * _uStride - _uSampStep);
			_uSampStep = 0;
			m2 = QDistri(QDis, cQ + 1);
			m1 = cQ * QDis[0];
			m1 = (m1 / nN - 1) / (cQ - 1);
			_fQuan << Arg_io[AI_HAMI] << endl;// '\t' << m2 << endl;
			_uCurSample++;
			bContinue = (_uCurSample < _uSamples);
			if (!bContinue) {
				sprintf(
						filename,
						"togcep%d_%d_%d~%d~%.15lg_%.15lg_%.15lg_%ld_%ldLastC.txt",
						nDim0, nDim1, int(cQ), nInitconf, dBeta, dAlpha, dUn,
						_uSamples, _uStride);
				DumpConf(filename);
			}
			if (_uCurSample % (_uSamples / 12) == 0) {
				outAcc = true;
			}

		}

		MetroplisTrialGCE(_uOneTimeMoves);
		_uSampStep += _uOneTimeMoves;

		//XArg_out[MCS_XAI_XKEY] = dUn;
		//XArg_out[MCS_XAI_HAMI] = Arg_io[AI_HAMI];//instant E

		return bContinue;
	}
	virtual bool TryExchange(TemperOMP_GCEPotts * XNew, double r)//random states
	{
		double x0 = -dAlpha * (dUn - XNew->dUn) * (Arg_io[AI_HAMI]
				- XNew->Arg_io[AI_HAMI]);
		double x = exp(x0);
		bool s = (r < x);

		if (s) {
			swap(States_p, XNew->States_p);
			swap(nCurHami, XNew->nCurHami);
			uAccepted++;//only add to the smalle one.
		}
		uChecked++;
		return s;
	}
};

void para_temper_SEO(TemperOMP_GCEPotts*mycomms, int nm, int doflag) {
	/*Stochastic even/odd algorithm(SEO),Chemical physics Letters 478(2009)80-84
	 * Potts 16*16,Q=10,\beta_0=1.420,E/N=-17568~-0.8687,N_rep=34, random next neighbor(RNN) 13.69 times of SEO*/

	unsigned long uOneTimeMoves = 1; //< moves interval between replica exchanges,"_uOneTimeMoves" is a random number
	int nRank, lk;
	bool bEvenIndexTemper;
	bool continued;
	Random *gran = new Random;
	double _dExchangedRate = nm - 1;
	_dExchangedRate /= mycomms[0].nN;//< set exchange rate 1%nProSize
	cout << "SEO_ExchangedRate:" << _dExchangedRate << '\t' << mycomms[0].nN
			- 1 << endl;

#pragma omp parallel if(nm > 1) shared(mycomms,nm,  gran,uOneTimeMoves,bEvenIndexTemper),private(nRank,continued,lk)
	{
		nRank = omp_get_thread_num();

		mycomms[nRank].Gauge();
		mycomms[nRank].MetroplisTrialGCE(-doflag * mycomms[nRank].nN);
		mycomms[nRank]._OpenOutPut(std::ios_base::out);
		continued = true;
		uOneTimeMoves = 1;
		while (continued) {

			continued = mycomms[nRank].RunMetroplisTrialGCE(uOneTimeMoves);
#pragma omp single
			{
				bEvenIndexTemper = gran->Real() > 0.5;
				//uNewIndexTemper = uOldIndexTemper + 1;
			}
			/* a barrier is automatically inserted here */
#pragma omp single
			{
				if (bEvenIndexTemper) {
					for (lk = 0; lk + 1 < nm; lk += 2)
						mycomms[lk].TryExchange(&mycomms[lk + 1], gran->Real());
				} else {
					for (lk = 1; lk + 1 < nm; lk += 2)
						mycomms[lk].TryExchange(&mycomms[lk + 1], gran->Real());
				}
				//Geometric distribution P(X=k)=(1-p)^{k-1}p, n trial: E(n)=1/p, var(n)=(1-p)/p^2
				//Uniform E(n)=(b+a)/2, var(n)=(b-a)^2/12;
				uOneTimeMoves = gran->Number(1, mycomms[0].nN - 1);
				//cout << uOldIndexTemper << '\t' << uOneTimeMoves << endl;
				if (mycomms[0].outAcc) {
					for (lk = 0; lk < nm - 1; lk++) {
						cout << mycomms[lk].dUn << '\t'
								<< mycomms[lk].uAccepted << '\t'
								<< double(mycomms[lk].uAccepted)
										/ double(mycomms[lk].uChecked) << endl;
						mycomms[lk].uAccepted = 0;
						mycomms[lk].uChecked = 0;
						mycomms[lk].outAcc = false;
					}
					cout << endl;
				}
			}
			//barrier
		}//while
		mycomms[nRank]._CloseOutPut();
	}
	delete gran;
}
int main(int argc, char *argv[]) {

	///check and set num_threads
	int nm = atoi(argv[11]);
	if (nm < 2) {
		cout << "please set 2<=np:" << nm << endl;
		return -1;
	} else {
		omp_set_num_threads(nm);
		int setok = false;
#pragma omp parallel shared(nm,setok)
		{
			setok = (nm == omp_get_num_threads());
			if (!setok) {
				cout << "omp_set_num_threads:error nm=" << nm << ':'
						<< omp_get_num_threads() << endl;
			}
		}
		if (!setok) {
			return -1;
		}
	}
	///task
	if (argc == 13) {
		if (nm != atoi(argv[11])) {///< we need a certain 'np'(number of processes),depending on the size of t[]
			cout << "please set np=" << nm << endl;
			return -1;
		}
		clock_t __start = clock();
		double *un0 = new double[nm];
		double unl = atof(argv[9]);
		double unh = atof(argv[10]);
		double beta0 = atof(argv[7]);
		double alpha = atof(argv[8]);
		unsigned long samples = atoi(argv[5]);
		unsigned long stride = atoi(argv[6]);
		int doflag = atoi(argv[12]);

		for (int i = 0; i < nm - 1; i++)//Linspace
			un0[i] = unl + i * (unh - unl) / (floor(nm) - 1);
		un0[nm - 1] = unh;
		int dim0 = atoi(argv[1]), dim1 = atoi(argv[2]), q = atoi(argv[3]),
				initconf = atoi(argv[4]);

		TemperOMP_GCEPotts *mycomms = new TemperOMP_GCEPotts[nm];

		if (doflag > 0) {
			int k;
#pragma omp parallel for shared(nm)
			for (k = 0; k < nm; k++) {
				mycomms[k].RunStatQuan(dim0, dim1, q, initconf, beta0, alpha,
						un0[k], samples, stride, doflag, k);// beta0, alpha,un0,samples, stride
			}
			cout << "Elapsed: " << ((double) (clock() - __start))
					/ CLOCKS_PER_SEC << std::endl;
		} else {
			for (int k = 0; k < nm; k++) {
				mycomms[k].Initiate(dim0, dim1, q, initconf, beta0, alpha,
						un0[k], samples, stride, doflag, k);
			}
			para_temper_SEO(mycomms, nm, doflag);
		}

		delete[] un0;
		delete[] mycomms;
		cout << "Elapsed: " << ((double) (clock() - __start)) / CLOCKS_PER_SEC
				<< std::endl;
	} else {
		cout
				<< "<TemperOMP_GCEPotts> dim0 dim0 Q initConf samples stride beta0,alpha,uNlow,uNhigh,nm,doflag"
				<< endl;
	}
	return 1;
}

