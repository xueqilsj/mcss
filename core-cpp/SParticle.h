/*
 * SParticle.h
 *  1D,U(x)
 *  Created on: Dec 12, 2008
 */

#ifndef _SParticle_H
#define _SParticle_H

#ifdef RANDOM_MM
#include "Random_mm.h"
#elif defined(RANDOM_CL)
#include "Random_cl.h"
#else
#include "Random.h"
#endif

class SParticle {
protected:

public:
	enum {
		AI_TEMPER = 0, AI_HAMI, AI_DELT_HAMI,
		///del.
		AI_MAG,
		AI_DELT_MAG,
		AI_RANDOM,
		AI_EXP_DELT_HAMI,
		AI_NUM
	//the last one define the element number of this enum
	};//derived class can extend this enum
	int nAccepted;
	double dLastTrial;
	double dX;// Point at
	double dOffset;//

	unsigned long long uCurMove;

	double Arg_io[AI_NUM];
	Random ran;//< random number generator
	SParticle();
	virtual ~SParticle();

	void SetArg(int site, double df);
	double GetArg(int site);

	double MP_XL;
	double MP_XH;
	void Gauge_MP();//multi point

	double PAK_XL;
	double PAK_XH;
	void Gauge_Pak();

	//double Delt();
	//========== NVT ===========//
	void InitParameters(double initPt, double off, double T);
	void MetroplisTrial_MP(int Moves);
	void MetroplisTrial_Pak(int Moves);

};

#endif
