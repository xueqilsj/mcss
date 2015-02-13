#include <omp.h>
#include "Potts.h"
#include "nvector.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

class TemperOMP_RCEPotts: public Potts {
protected:
	ofstream _fQuan;
	unsigned long _uSamples, _uCurSample, _uStride, _uSampStep;
	int nRank;
	///////////
public:
	char filename[120];
	int nInitconf;
	double UN0, Alpha;
	unsigned long uAccepted;
	unsigned long uChecked;
	bool outAcc;
	TemperOMP_RCEPotts() {
	}
	virtual ~TemperOMP_RCEPotts() {
	}
	void _OpenOutPut(double un0, std::ios::openmode mode) {
		sprintf(filename,
				"tpt%d_%d_%d~%d~%lf_%.14lf_%lf_%ld_%ldQuanFull.txt", nDim0,
				nDim1, int(cQ), nInitconf, Arg_io[AI_TEMPER], un0, Alpha,
				_uSamples, _uStride);
		cout << filename << endl;
		_fQuan.open(filename, mode);
	}
	void _CloseOutPut(double un0) {
		_fQuan.close();
	}
	///////////
	int Setting(double un0, double dBeta0, double alpha, unsigned long Samples,
			unsigned long Stride, int Rank) {

		if (Stride == 0) {
			cout << "please set _uStride nonzero" << endl;
			return -1;
		}
		_uSamples = Samples;
		_uStride = Stride;
		_uSampStep = 0;
		_uCurSample = 0;
		nRank = Rank;
		uAccepted = 0;
		uChecked = 0;
		outAcc = false;
		UN0 = un0;
		Arg_io[AI_TEMPER] = dBeta0;
		Alpha = alpha;

		time_t now = time(NULL);
		now /= 31;
		now *= 18653323 + nRank * 1231532129; //< Let every replica different random seed
		ran.Srand(now);
		return 1;
	}
	int LoadLastC(double un0, double dBeta0, double alpha,
			unsigned long Samples, unsigned long Stride) {

		if (Stride == 0) {
			cout << "please set _uStride nonzero" << endl;
			return -1;
		}
		sprintf(filename, "tpt%d_%d_%d~%d~%lf_%.14lf_%lf_%ld_%ldLastC.txt",
				nDim0, nDim1, int(cQ), nInitconf, dBeta0, un0, alpha, Samples,
				Stride);

		if (!LoadConf(filename)) {
			cerr << "#LoadConf Error:" << filename << endl;
			return -1;
		}
		return 1;
	}
	int RunStatQuan(double un0, double dBeta0, double alpha,
			unsigned long Samples, unsigned long Stride, unsigned long rehead,
			int nRank) {

		if (Stride == 0) {
			cout << "please set _uStride nonzero" << endl;
			return -1;
		}
		_uSamples = Samples;
		_uStride = Stride;

		UN0 = un0;
		Arg_io[AI_TEMPER] = dBeta0;
		Alpha = alpha;
		sprintf(filename,
				"tpt%d_%d_%d~%d~%lf_%.14lf_%lf_%ld_%ldQuanFull.txt", nDim0,
				nDim1, int(cQ), nInitconf, Arg_io[AI_TEMPER], UN0, Alpha,
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
		double mean, std, mean1, std1;

		unsigned int neff = Esamples.mean_std(&mean, &std, false, uTdrop, uTau);
		double Beta = Arg_io[AI_TEMPER] + Alpha * (mean / nN - UN0);
		Esamples.mean_std(&mean1, &std1, false, uTdrop, 1);

		sprintf(filename, "tpt%d_%d_%d~%d~%lf_%02d_%lf_%ld_%ldStat.txt",
				nDim0, nDim1, int(cQ), nInitconf, Arg_io[AI_TEMPER], nRank,
				Alpha, _uSamples, _uStride);
		ofstream _fStatOut(filename);
		_fStatOut << "#nInit\t" << setw(9) << "UN0"
				<< "\tB0\tAlpha\tWin\tN\tTeq\tnLags\tTau\tNeff\tmean\tstd\tBeta\tmean1\tstd1"
				<< endl << nInitconf << '\t' << setw(9) << setprecision(16)
				<< UN0 << '\t' << Arg_io[AI_TEMPER] << '\t' << Alpha << '\t'
				<< win << '\t' << Esamples.len << '\t' << uTdrop << '\t'
				<< uLags << '\t' << uTau << '\t' << neff << '\t' << mean
				<< '\t' << std << '\t' << Beta << '\t' << mean1 << '\t' << std1
				<< endl;
		_fStatOut.close();
		return 1;
	}
	virtual bool _NVTDo(unsigned long _uOneTimeMoves) {
		bool bContinue = true;
		unsigned int QDis[cQ + 1];
		double m1, m2;

		while ((_uSampStep + _uOneTimeMoves) >= nN * _uStride) {//sampling
			MetroplisTrialRCE(nN * _uStride - _uSampStep, UN0,
					Arg_io[AI_TEMPER], Alpha);
			_uOneTimeMoves -= (nN * _uStride - _uSampStep);
			_uSampStep = 0;
			m2 = QDistri(QDis, cQ + 1);
			m1 = cQ * QDis[0];
			m1 = (m1 / nN - 1) / (cQ - 1);
			_fQuan << setprecision(16) << Arg_io[AI_HAMI] << '\t'
					<< setprecision(6) << m1 << endl;// '\t' << m2 << endl;
			_uCurSample++;
			bContinue = (_uCurSample < _uSamples);
			if (!bContinue) {
				sprintf(filename,
						"tpt%d_%d_%d~%d~%lf_%.14lf_%lf_%ld_%ldLastC.txt",
						nDim0, nDim1, int(cQ), nInitconf, Arg_io[AI_TEMPER],
						UN0, Alpha, _uSamples, _uStride);
				DumpConf(filename);
			}
			if (_uCurSample % (_uSamples / 20) == 0) {
				outAcc = true;
			}

		}

		MetroplisTrialRCE(_uOneTimeMoves, UN0, Arg_io[AI_TEMPER], Alpha);
		_uSampStep += _uOneTimeMoves;

		//XArg_out[MCS_XAI_XKEY] = UN0;
		//XArg_out[MCS_XAI_HAMI] = Arg_io[AI_HAMI];//instant E

		return bContinue;
	}
	virtual bool TryExchange(TemperOMP_RCEPotts * XNew, double r)//random states
	{
		double x0 = -Alpha * (UN0 - XNew->UN0) * (Arg_io[AI_HAMI]
				- XNew->Arg_io[AI_HAMI]);
		double x = exp(x0);
		bool s = (r < x);

		if (s) {
			swap(States_p, XNew->States_p);
			swap(nCurHami, XNew->nCurHami);
			uAccepted++;
		}
		uChecked++;
		return s;
	}
};

void para_temper_SEO(TemperOMP_RCEPotts*mycomms, int nm, double *un0,
		double beta0, double alpha, int doflag) {
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

#pragma omp parallel if(nm > 1) shared(mycomms,nm, un0, beta0, alpha, gran,uOneTimeMoves,bEvenIndexTemper),private(nRank,continued,lk)
	{
		nRank = omp_get_thread_num();

		mycomms[nRank].Gauge();
		mycomms[nRank].MetroplisTrialRCE(-doflag * mycomms[nRank].nN,
				un0[nRank], beta0, alpha);
		mycomms[nRank]._OpenOutPut(un0[nRank], std::ios_base::out);
		continued = true;
		uOneTimeMoves = 1;
		while (continued) {

			continued = mycomms[nRank]._NVTDo(uOneTimeMoves);
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
				if (mycomms[nRank].outAcc) {
					for (lk = 0; lk < nm; lk++) {
						cout << mycomms[lk].UN0 << '\t'
								<< mycomms[lk].uAccepted << '\t'
								<< double(mycomms[lk].uAccepted)
										/ double(mycomms[lk].uChecked) << endl;
						mycomms[lk].outAcc = false;
					}
					cout << endl;
				}
			}
			//barrier
		}//while
		mycomms[nRank]._CloseOutPut(un0[nRank]);
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
				cout << "omp_set_num_threads:error nm=" << nm << endl;
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

		for (int i = 0; i < nm - 1; i++)
			un0[i] = unl + i * (unh - unl) / (floor(nm) - 1);
		un0[nm - 1] = unh;
		int dim0 = atoi(argv[1]), dim1 = atoi(argv[2]), q = atoi(argv[3]),
				iniconf = atoi(argv[4]);
		char *States_ps = new char[nm * (dim0 * dim1) + 1];
		TemperOMP_RCEPotts *mycomms = new TemperOMP_RCEPotts[nm];

		for (int k = 0; k < nm; k++) {
			mycomms[k].InitParameters(q);//< Dimension=2,Q=10
			mycomms[k].SetConf(States_ps + k * (dim0 * dim1), dim0, dim1);
			mycomms[k].nInitconf = iniconf;
			if (doflag == 0) {
				if (mycomms[k].LoadLastC(un0[k], beta0, alpha, samples, stride)
						== -1)
					return -1;
			} else if (doflag < 0) {
				mycomms[k].InitConf(mycomms[k].nInitconf);
			}
			mycomms[k].Setting(un0[k], beta0, alpha, samples, stride, k);//un0, beta0, alpha,samples, stride

		}
		if (doflag > 0) {
			int k;
#pragma omp parallel for shared(nm)
			for (k = 0; k < nm; k++) {
				mycomms[k].RunStatQuan(un0[k], beta0, alpha, samples, stride,
						doflag, k);//un0, beta0, alpha,samples, stride
			}
			cout << "Elapsed: " << ((double) (clock() - __start))
					/ CLOCKS_PER_SEC << std::endl;
		} else {
			para_temper_SEO(mycomms, nm, un0, beta0, alpha, doflag);
		}
		delete[] un0;
		delete[] States_ps;
		delete[] mycomms;
		cout << "Elapsed: " << ((double) (clock() - __start)) / CLOCKS_PER_SEC
				<< std::endl;
	} else {
		cout
				<< "<TemperOMP_RCEPotts> dim0 dim0 Q initConf samples stride beta0,alpha,uNlow,uNhigh,nm,doflag"
				<< endl;
	}
	return 1;
}

