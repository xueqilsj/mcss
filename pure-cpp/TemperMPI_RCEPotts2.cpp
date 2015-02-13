/*
 ============================================================================
 Name        : TemperMPI_RCEPotts.cpp
 Description : Sampling in Potts model with rotational canonical Ensemble:beta(E)=beta0+alpha*(E-E0)/N.
 ============================================================================
 */

#include "Potts.h"
#include "MPITemper.h"
#include "LogFile.h"
#include "nvector.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

class TemperMPI_RCEPotts: public Potts, public MPITemper {
protected:
	LogFile _fQuan;
	double _dAcceptedRate, UN0, Alpha;
	unsigned long _uSamples, _uCurSample, _uStride;
	long _uSampStep;
	int nInitconf;
	///////////
public:
	char filename[120];
	TemperMPI_RCEPotts(int dim0, int dim1, char Q, int argc, char * argv[],
			int initconf = 0) :
		MPITemper(argc, argv) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		nInitconf = initconf;
		InitConf(nInitconf);

	}
	virtual ~TemperMPI_RCEPotts() {
		delete[] States_p;
	}

	void _OpenOutPut(double un0, std::ios::openmode mode) {
		sprintf(filename,
				"trcep%d_%d_%d~%d~%lf_%.14lf_%lf_%ld_%ldQuanFull.txt", nDim0,
				nDim1, int(cQ), nInitconf, Arg_io[AI_TEMPER], un0, Alpha,
				_uSamples, _uStride);
		_fQuan.Open(filename, mode);
	}
	void _CloseOutPut(double un0) {
		_fQuan.Close();
	}
	///////////
	int Apportion(int nNumOfUNTable, double dUNTable[], double dBeta0,
			double alpha, unsigned long Samples, unsigned long Stride) {
		if (nProSize != nNumOfUNTable) {///< we need a certain 'np'(number of processes),depending on the size of t[]
			LOG_ERROR2(cout, "please set np=", nNumOfUNTable);
			return -1;
		}
		if (Stride == 0) {
			LOG_ERROR(cout, "please set _uStride nonzero");
			return -1;
		}
		_uSamples = Samples;
		_uStride = Stride;
		_uSampStep = 0;
		_uCurSample = 0;
		_dExchangedRate = nProSize * 1.0 / 1000.0;//< set exchange rate 1%nProSize
		_uOneTimeMoves = 1;
		Gauge();
		MPITemper::_ParaInit(nNumOfUNTable, dUNTable);
		UN0 = _TemperTable_p[nRank].dXKey;
		Arg_io[AI_TEMPER] = dBeta0;
		Alpha = alpha;

		_OpenOutPut(UN0, std::ios_base::out);

		time_t now = time(NULL);
		now /= 31;
		now *= 18653323 + nRank * 1231532129; //< Let every replica different random seed
		ran.Srand(now);

		int res = MPITemper::Run();

		if (nRank == 0) {
			for (int i = 0; i < nProSize - 1; i++)
				cout << _TemperTable_p[i].dXKey << '\t'
						<< _TemperTable_p[i].uAccepted << '\t'
						<< double(_TemperTable_p[i].uAccepted)
								/ double(_TemperTable_p[i].uChecked) << endl;
			cout << endl;
		}
		_CloseOutPut(UN0);
		return res;
	}
	int LoadLastC(int nNumOfUNTable, double dUNTable[], double dBeta0,
			double alpha, unsigned long Samples, unsigned long Stride) {
		if (nProSize != nNumOfUNTable) {///< we need a certain 'np'(number of processes),depending on the size of t[]
			LOG_ERROR2(cout, "please set np=", nNumOfUNTable);
			return -1;
		}
		if (Stride == 0) {
			LOG_ERROR(cout, "please set _uStride nonzero");
			return -1;
		}
		sprintf(filename, "trcep%d_%d_%d~%d~%lf_%.14lf_%lf_%ld_%ldLastC.txt",
				nDim0, nDim1, int(cQ), nInitconf, dBeta0, dUNTable[nRank],
				alpha, Samples, Stride);

		if (!LoadConf(filename)) {
			cerr << "#LoadConf Error:" << filename << endl;
			return -1;
		}
		return 1;
	}
	int RunStatQuan(int nNumOfUNTable, double dUNTable[], double dBeta0,
			double alpha, unsigned long Samples, unsigned long Stride,
			unsigned long rehead) {
		if (nProSize != nNumOfUNTable) {///< we need a certain 'np'(number of processes),depending on the size of t[]
			LOG_ERROR2(cout, "please set np=", nNumOfUNTable);
			return -1;
		}
		if (Stride == 0) {
			LOG_ERROR(cout, "please set _uStride nonzero");
			return -1;
		}
		_uSamples = Samples;
		_uStride = Stride;

		MPITemper::_ParaInit(nNumOfUNTable, dUNTable);
		UN0 = _TemperTable_p[nRank].dXKey;
		Arg_io[AI_TEMPER] = dBeta0;
		Alpha = alpha;
		sprintf(filename,
				"trcep%d_%d_%d~%d~%lf_%.14lf_%lf_%ld_%ldQuanFull.txt", nDim0,
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

		sprintf(filename, "trcep%d_%d_%d~%d~%lf_%02d_%lf_%ld_%ldStat.txt",
				nDim0, nDim1, int(cQ), nInitconf, Arg_io[AI_TEMPER], nRank,
				Alpha, _uSamples, _uStride);
		LogFile _fStatOut;
		_fStatOut.Open(filename, std::ios_base::out);
		LOG_STREAM(_fStatOut) << "#nInit\t" << setw(9) << "UN0"
				<< "\tB0\tAlpha\tWin\tN\tTeq\tnLags\tTau\tNeff\tmean\tstd\tBeta\tmean1\tstd1"
				<< endl << nInitconf << '\t' << setw(9) << setprecision(16)
				<< UN0 << '\t' << Arg_io[AI_TEMPER] << '\t' << Alpha << '\t'
				<< win << '\t' << Esamples.len << '\t' << uTdrop << '\t'
				<< uLags << '\t' << uTau << '\t' << neff << '\t' << mean
				<< '\t' << std << '\t' << Beta << '\t' << mean1 << '\t' << std1
				<< endl;
		_fStatOut.Close();
		return 0;
	}
	virtual bool _NVTDo(double XArg_out[]) {
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
			LOG_STREAM(_fQuan) << setprecision(16) << Arg_io[AI_HAMI] << '\t'
					<< setprecision(6) << m1 << '\t' << m2 << endl;
			_uCurSample++;
			bContinue = (_uCurSample < _uSamples);
			if (!bContinue) {
				sprintf(filename,
						"trcep%d_%d_%d~%d~%lf_%.14lf_%lf_%ld_%ldLastC.txt",
						nDim0, nDim1, int(cQ), nInitconf, Arg_io[AI_TEMPER],
						UN0, Alpha, _uSamples, _uStride);
				DumpConf(filename);
			}
			if (nRank == 0 && _uCurSample % _uSamples / 5 == 0) {
				for (int i = 0; i < nProSize - 1; i++)
					cout << _TemperTable_p[i].dXKey << '#'
							<< _TemperTable_p[i].uAccepted << ':'
							<< double(_TemperTable_p[i].uAccepted)
									/ double(_TemperTable_p[i].uChecked) << ' ';
				cout << endl;
			}
		}

		MetroplisTrialRCE(_uOneTimeMoves, UN0, Arg_io[AI_TEMPER], Alpha);
		_uSampStep += _uOneTimeMoves;

		XArg_out[MCS_XAI_XKEY] = UN0;
		XArg_out[MCS_XAI_HAMI] = Arg_io[AI_HAMI];//instant E
		XArg_out[MCS_XAI_ACCEPTED_RATE] = _dAcceptedRate;
		XArg_out[MCS_XAI_CUR_TINDEX] = _uCurIndexTemper;
		// TODO: At present,don't transfer _dMagneticAll, _dMagneticAll2
		return bContinue;
	}
	virtual void _UpdateExchange(double XOtherArg_in[]) {
		UN0 = XOtherArg_in[MCS_XAI_XKEY];
		Arg_io[AI_HAMI] = XOtherArg_in[MCS_XAI_HAMI];
		_dAcceptedRate = XOtherArg_in[MCS_XAI_ACCEPTED_RATE];

	}
	virtual bool _CheckExchange(double XOldArg_in[], double XNewArg_in[])//random states
	{
		double x0 = -Alpha * (XOldArg_in[MCS_XAI_XKEY]
				- XNewArg_in[MCS_XAI_XKEY]) * (XOldArg_in[MCS_XAI_HAMI]
				- XNewArg_in[MCS_XAI_HAMI]);
		double x = exp(x0);
		double r = _GetRandomReal(0, 1.0);
		bool s = (r < x);
		return s;
	}
	virtual inline double _GetRandomReal(double low, double high) {
		return ran.Real(low, high);
	}
	virtual inline int _GetRandomNumber(int low, int high) {
		return ran.Number(low, high);
	}

};

int main(int argc, char *argv[]) {

	if (argc == 13) {
		int nm = atoi(argv[11]);
		double *un0 = new double[nm];
		double unl = atof(argv[9]);
		double unh = atof(argv[10]);
		for (int i = 0; i < nm - 1; i++)
			un0[i] = unl + i * (unh - unl) / (floor(nm) - 1);
		un0[nm - 1] = unh;

		TemperMPI_RCEPotts mycomm(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]),
				argc, argv, atoi(argv[4]));//< Dimension=2,Q=10
		double beta0 = atof(argv[7]);
		double alpha = atof(argv[8]);
		unsigned long samples = atoi(argv[5]);
		unsigned long stride = atoi(argv[6]);
		int doflag = atoi(argv[12]);

		if (doflag < 0) {
			mycomm.Gauge();
			mycomm.MetroplisTrialRCE(-doflag * mycomm.nN, un0[mycomm.nRank],
					beta0, alpha);
			mycomm.Apportion(nm, un0, beta0, alpha, samples, stride);//un0, beta0, alpha,samples, stride
		} else if (doflag == 0) {
			if (mycomm.LoadLastC(nm, un0, beta0, alpha, samples, stride) != -1) {
				mycomm.Apportion(nm, un0, beta0, alpha, samples, stride);//un0, beta0, alpha,samples, stride
			}
		} else {
			mycomm.RunStatQuan(nm, un0, beta0, alpha, samples, stride, doflag);//un0, beta0, alpha,samples, stride
		}
	} else {
		cout
				<< "<TemperMPI_RCEPotts> dim0 dim0 Q initConf samples stride beta0,alpha,uNlow,uNhigh,nm,doflag"
				<< endl;
	}
}
