/*
 * Linear.h
 *
 *  Created on: Dec 12, 2008
 */

#ifndef _Linear_H
#define _Linear_H

#include "Random.h"

class Linear {
protected:

public:
	enum {
		AI_TEMPER = 0,
		AI_HAMI,
		AI_DELT_HAMI,
		AI_MAG,
		AI_DELT_MAG,
		AI_RANDOM,
		AI_EXP_DELT_HAMI,
		AI_NUM
	//the last one define the element number of this enum
	};//derived class can extend this enum
	int nAccepted;
	double dLastTrial[1];
	double dDim[2];//< (l1,h1,l2,h2,..lD,HD)
	double dPoint[1];//< Point=(x1,x2,x3,...)
	double dOffset[1];//< Point=(x1,x2,x3,...)

	unsigned long long uCurMove;

	double Arg_io[AI_NUM];
	Random ran;//< random number generator
	Linear();
	virtual ~Linear();

	void SetArg(int site, double df);
	double GetArg(int site);
	double GetPoint(int in);

	void Gauge();
	double Delt();
	//========== NVT ===========//
	void InitParameters(double initp, double off, double T);
	int Accepted();//NVT
	void MetroplisTrial(int Moves);

};

#endif
