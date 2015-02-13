/*
 * MPITemper.cpp
 *
 *  Created on: Jan 7, 2009
 */

#include "MPITemper.h"
#include "LogFile.h"
#include <cmath> //for exp()
using namespace std;

MPITemper::MPITemper(int argc, char * argv[]) :
	_uOneTimeMoves(_uRefArg_a[0]), _uOldIndexTemper(_uRefArg_a[1]),
			_uNewIndexTemper(_uRefArg_a[2]), _uLastTimeCheckOK(_uRefArg_a[3]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProSize);
	_TemperTable_p = NULL;
}

MPITemper::~MPITemper() {
	MPI_Finalize();
	if (_TemperTable_p != NULL)
		delete[] _TemperTable_p;
}
int MPITemper::Run(unsigned long XAINum) {
	_uMcsXAINum = XAINum;
	//_ParaInit(nNumOfTemperatureTable, dTemperatureTable);

	double starttime = MPI_Wtime();
	double XArg_iot[_uMcsXAINum];//< variables for replica exchange: T,E,(_dEnergyAll)/AVERAGE_Energy
	unsigned long uLastOldIndexTemper_t = _uOldIndexTemper;
	unsigned long uLastNewIndexTemper_t = _uNewIndexTemper;
	while (true) {
		_BCastAndReceive(uLastOldIndexTemper_t, uLastNewIndexTemper_t);
		if (!_NVTDo(XArg_iot))
			break;
		uLastOldIndexTemper_t = _uOldIndexTemper;
		uLastNewIndexTemper_t = _uNewIndexTemper;

		//_uCurIndexTemper may be changed above!
		if (_uCurIndexTemper == _uOldIndexTemper) {
			_OldDo(XArg_iot);
		} else if (_uCurIndexTemper == _uNewIndexTemper) {
			_NewDo(XArg_iot);
		}
	}

	cout << "rank: " << nRank << " exit, eclapsed=" << MPI_Wtime() - starttime
			<< std::endl;
	return nRank;
}
void MPITemper::_ParaInit(int nNumOfXKeyTable, double dXkeyTable[]) {
	_uLastTimeCheckOK = 0;
	_uCurIndexTemper = nRank;
	_uOldIndexTemper = 0;
	_uNewIndexTemper = 1;

	_TemperTable_p = new struct Temper_st[nNumOfXKeyTable];
	for (int i = 0; i < nNumOfXKeyTable; i++) {
		_TemperTable_p[i].dXKey = dXkeyTable[i];
		_TemperTable_p[i].uProcess = i;
		_TemperTable_p[i].uAccepted = 0;
		_TemperTable_p[i].uChecked = 0;

	}

}
void MPITemper::_BCastAndReceive(const unsigned long &lastOldIndexTemper,
		const unsigned long &lastNewIndexTemper) {
	// _nMoves, oldIndexT, newIndexT,	lastChOK,
	MPI_Bcast(_uRefArg_a, 4, MPI_UNSIGNED_LONG,
			_TemperTable_p[lastOldIndexTemper].uProcess, MPI_COMM_WORLD);//blocking receive
	///_uRefArg_a keep unchanged below
	if (_uLastTimeCheckOK != 0) {//At last time, checking is ok
		swap(_TemperTable_p[lastOldIndexTemper].uProcess,
				_TemperTable_p[lastNewIndexTemper].uProcess);

		if (lastOldIndexTemper < lastNewIndexTemper)
			_TemperTable_p[lastOldIndexTemper].uAccepted++;
		else
			_TemperTable_p[lastNewIndexTemper].uAccepted++;
	}

	if (lastOldIndexTemper < lastNewIndexTemper)
		_TemperTable_p[lastOldIndexTemper].uChecked++;
	else
		_TemperTable_p[lastNewIndexTemper].uChecked++;

}

void MPITemper::_OldDo(double XOldArg_io[])//<from _NVT()
{
	unsigned long uNewIndexTemper_t = _uNewIndexTemper;//<due to _uNewIndexTemper will be modified below!
	MPI_Send(XOldArg_io, _uMcsXAINum, MPI_DOUBLE,
			_TemperTable_p[uNewIndexTemper_t].uProcess, 2, MPI_COMM_WORLD);
	_CloseOutPut(_TemperTable_p[_uCurIndexTemper].dXKey);
	_uOneTimeMoves = 0;
	while (_dExchangedRate < _GetRandomReal(0, 1))		_uOneTimeMoves++;//binomial distribution
	//_uOneTimeMoves = _GetRandomNumber(0, (1.0 - _dExchangedRate) * 200);
	_uOldIndexTemper = _GetRandomNumber(0, nProSize - 1);

	if (_uOldIndexTemper == static_cast<unsigned long> (nProSize - 1)) {
		_uNewIndexTemper = _uOldIndexTemper - 1;
	} else if (_uOldIndexTemper == 0) {
		_uNewIndexTemper = _uOldIndexTemper + 1;
	} else {
		_uNewIndexTemper = (_GetRandomNumber(0, 1) == 0) ? _uOldIndexTemper - 1
				: _uOldIndexTemper + 1;
	}

	MPI_Recv(XOldArg_io, _uMcsXAINum, MPI_DOUBLE,
			_TemperTable_p[uNewIndexTemper_t].uProcess, 3, MPI_COMM_WORLD,
			&_Status);//blocking until received T,E
	_uLastTimeCheckOK = (_uCurIndexTemper
			!= static_cast<unsigned long> (XOldArg_io[MCS_XAI_CUR_TINDEX]));
	if (_uLastTimeCheckOK == 1) {
		_UpdateExchange(XOldArg_io);
		_uCurIndexTemper = uNewIndexTemper_t;
	}
	_OpenOutPut(_TemperTable_p[_uCurIndexTemper].dXKey, std::ios_base::out
			| std::ios_base::app);

}

void MPITemper::_NewDo(double XNewArg_io[])//<from _NVT()
{
	double XOldArg_iot[_uMcsXAINum];
	bool bEXchang;
	MPI_Recv(XOldArg_iot, _uMcsXAINum, MPI_DOUBLE,
			_TemperTable_p[_uOldIndexTemper].uProcess, 2, MPI_COMM_WORLD,
			&_Status);//blocking until received T,E
	bEXchang = _CheckExchange(XOldArg_iot, XNewArg_io);
	if (bEXchang) {
		_CloseOutPut(_TemperTable_p[_uCurIndexTemper].dXKey);
		_UpdateExchange(XOldArg_iot);
		for (unsigned int n = 0; n < _uMcsXAINum; n++) {
			swap(XOldArg_iot[n], XNewArg_io[n]);
		}
		_uCurIndexTemper = _uOldIndexTemper;
	}
	MPI_Send(XOldArg_iot, _uMcsXAINum, MPI_DOUBLE,
			_TemperTable_p[_uOldIndexTemper].uProcess, 3, MPI_COMM_WORLD);//T,E
	if (bEXchang)
		_OpenOutPut(_TemperTable_p[_uCurIndexTemper].dXKey, std::ios_base::out
				| std::ios_base::app);
}

bool MPITemper::_CheckExchange(double XOldArg_in[], double XNewArg_in[])//random states
{
	double x0 = (1.0 / XOldArg_in[MCS_XAI_XKEY] - 1.0
			/ XNewArg_in[MCS_XAI_XKEY]) * (XOldArg_in[MCS_XAI_HAMI]
			- XNewArg_in[MCS_XAI_HAMI]);
	double x = exp(x0);
	double r = _GetRandomReal(0, 1.0);
	bool s = (r < x);

	return s;
}

