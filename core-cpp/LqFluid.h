/*
 * LqFluid.h
 *
 *  Created on: Dec 12, 2008
 *
 *   0---N-->
 *   ----x,dim0--->
 *   |
 *   |
 *   y,      (0,0)
 *  dim0
 *   |
 *   v
 *
 * 0,1,...,4,nDim1-1
 * nDim1,...2*nDim1-1
 * ...
 * (nDim0-1)*nDim1,...nDim0*nDim1-1
 *
 * Energy Levels(2N-3):
 * -2N,_,_,_,-2N-4,_,-2N-5,2N-6,..,0
 * 0,...4,,6
 *
 * first-order transition: Q>4, deltaE:-4,-3,-2,-1,1,0,1,2,3,4
 */

#ifndef _LqFluid_H
#define _LqFluid_H

#ifdef RANDOM_MM
#include "Random_mm.h"
#elif defined(RANDOM_CL)
#include "Random_cl.h"
#else
#include "Random.h"
#endif

class LqFluid {
protected:

public:
	double _dDeltaExp[9];
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
	double Arg_io[AI_NUM];
	double Evdw, delEvdw;
	unsigned int nDim0, nDim1, nN, nCurTrial;// Dimensions and size
	double pbcx, pbcy;//Periodic boundary condition
	double *Xp, *Yp;//, double *Zp;size=nDim0*nDim1
	double *FXp, *FYp;//, double *Zp;
	double CorVar[3];//x,y,z
	double FVar[3];//x,y,z :force
	double DEvdw[2];//old,new;

	double pairlistdist;// PAIRLISTDIST;
	double dDist, dMax, dt;//span,TIMEFACTOR
	double vdwA, vdwB, epsilon, sigma, cutoff, cut2, switchdist, switch2;
	double dAcceptRatio;
	double dBeta, dAlpha, dUn;//for Generalized canonical ensemble (dUn=U/N).
	int nInitconf;//starting configuration:
	//d<0 ground state;d==0 random state;d>0 anti-ground state,2(backup)
	//nN:size of model,nEn: the total number of energy levels
	unsigned long long uCurMove, uAccepted;
	int nEn, nCurHami, nDeltaE;//nEn:the number of total energy levels
	int nLowEBound, nHighEBound;//[nLowEBound, nHighEBound],for flat histogram:WL,Mul
	Random ran, ran2;//random number generator,automatic initiation.

	LqFluid();
	virtual ~LqFluid();

	void InitParameters(double ds, double dmax);
	void ResetLQ(double epsil, double sig, double swit, double cutf);
	double GetLqValue(double dist);
	int GetRDF(double *rdf, int len);
	void SetArg(int site, double df);
	double GetArg(int site);

	void SetConf(double *Corxs, double *Corys, double *Fxs, double *Fys,
			int dim0, int dim1);
	void InitConf(int d);//d>-2

	void Gauge();//calculate Hamiltonian of the current configuration.
	int Delt();//calculate delta Energy if the value of current spin changes to QVar.
	double pbc_length_ji2(int i, int j);
	void addF_ji(int i, int j, double f);

	void MetroplisSweepBeta(int sweeps);
	void MetroplisSweepGCE(int sweeps);

};
#endif
