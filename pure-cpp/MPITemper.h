/*
 * MPITemper.h
 *
 *  Created on: Jan 7, 2009
 */

#ifndef MPITEMPER_H_
#define MPITEMPER_H_

#include "mpi.h"
#include <iostream>
#include <list>

typedef unsigned long long UMove_type;//unsigned move
struct Temper_st {
	double dXKey;
	unsigned int uProcess;
	UMove_type uAccepted;
	UMove_type uChecked;
	int nExtDir;
};

class MPITemper {
protected:
	enum {//index for XArg[]
		MCS_XAI_XKEY = 0,
		MCS_XAI_HAMI,
		MCS_XAI_1,
		MCS_XAI_2,
		MCS_XAI_ACCEPTED_RATE,
		MCS_XAI_CUR_TINDEX,
		MCS_XAI_NUM
	};
	MPI_Status _Status;
	unsigned int _uMcsXAINum;
	struct Temper_st *_TemperTable_p;
	double _dExchangedRate; //<exchange rate for replica change
	unsigned long _uRefArg_a[4]; //< referenced by _uOneTimeMoves,_uOldIndexTemper,_uNewIndexTemper,_uLastTimeCheckOK
	unsigned long &_uOneTimeMoves; //< moves interval between replica exchanges,"_uOneTimeMoves" is a random number
	unsigned long &_uOldIndexTemper;
	unsigned long &_uNewIndexTemper;
	unsigned long &_uLastTimeCheckOK;
	unsigned long _uCurIndexTemper;
protected:
	virtual void _ParaInit(int nNumOfXkeyTable, double dXkeyTable[]);
	virtual void _BCastAndReceive(const unsigned long &lastOldIndexTemper,
			const unsigned long &lastNewIndexTemper);
	void _OldDo(double XArg_in[3]);
	void _NewDo(double XArg_in[3]);

	virtual bool _NVTDo(double XArg_out[])=0;
	virtual bool _CheckExchange(double XOldArg_in[], double XNewArg_in[]);
	virtual void _UpdateExchange(double XOtherArg_in[])=0;
	virtual double _GetRandomReal(double low, double high)=0;
	virtual int _GetRandomNumber(int low, int high)=0;
	virtual void _CloseOutPut(double xkey)=0;
	virtual void _OpenOutPut(double xkey, std::ios::openmode mode)=0;

public:
	int nRank, nProSize;
	MPITemper(int argc, char * argv[]);
	virtual int Run(unsigned long XAINum = MCS_XAI_NUM);

	virtual ~MPITemper();

};

#endif /* MPICOMM_H_ */
