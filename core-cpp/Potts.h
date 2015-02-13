/*
 * Potts.h
 *
 *  Created on: Dec 12, 2008
 *
 *   0---N-->
 * (0,0)----y,dim1--->
 *   |
 *   |
 *   x,
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

#ifndef _Potts_H
#define _Potts_H

#ifdef RANDOM_MM
#include "Random_mm.h"
#elif defined(RANDOM_CL)
#include "Random_cl.h"
#else
#include "Random.h"
#endif

class Potts {
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

	char QVar, cQ;//QVar in [1,cQ]
	unsigned int nDim0, nDim1, nN, nCurTrial;// Dimensions and size
	char * States_p;//configuration, size=nDim0*nDim1;
	double dBeta, dAlpha, dUn;//for Generalized canonical ensemble (dUn=U/N).
	int nInitconf;//starting configuration:
	//d<0 ground state;d==0 random state;d>0 anti-ground state,2(backup)
	//nN:size of model,nEn: the total number of energy levels
	unsigned long long uCurMove;
	int nEn, nCurHami, nDeltaE;//nEn:the number of total energy levels
	int nLowEBound, nHighEBound;//[nLowEBound, nHighEBound],for flat histogram:WL,Mul
	Random ran, ran2;//random number generator,automatic initiation.

	Potts();
	virtual ~Potts();

	void InitParameters(char q);
	void SetArg(int site, double df);
	double GetArg(int site);

	void SetConf(char * State_p, int dim0, int dim1);
	void InitConf(int d);//d<0 ground state;d==0 random state;d>0 anti-ground state
	bool LoadConf(const char * ConfName);
	bool FromConf(const char * ConfName);
	bool DumpConf(const char * ConfName);

	void Gauge();//calculate Hamiltonian of the current configuration.
	int Delt();//calculate delta Energy if the value of current spin changes to QVar.
	void DeltS(unsigned int CTrial);//TrialValue may be same as QVar
	//========== Order Parameter===========//
	double QDistri(unsigned int * QDis_p, int M);//M=Q+1,Nmax return Nvec2;
	//========== NVT Metroplis===========//
	void SetTemperature(double T);
	void MetroplisTrial(int Moves);//index exp(-dE)
	void MetroplisTrialBT(int Moves);//without index exp(-dE)
	//========== GCE Metroplis===========//
	void MetroplisTrialGCE(int Moves);
	//========== NVE Metroplis===========//
	double dRefU;//[-2N,0]==2N+1,Demon Energy
	bool SetRefU(double u);
	void MetroplisTrialNVE(int Moves);//Ray's NVE: \beta(E)=(N-2)/2(U-E),U is the total energy.
	int nEDemon;
	void MetroplisTrialCreutzNVE(int Moves);//Creutz
	//===========NVT Wolff Cluster=======//
	int * ClusterStack_p;
	double dPadd;//the probability of adding bonds

	void WolffCluster(int Moves);//
	int * ClusterSeq_p;//roll back if GCE cluster flipping refused.

	int csp;//the size of current cluster
	int nSegAccumCS;//accumulate size a part of cluster
	double dAccumCS;

	int nGCEWF_AvgC;//the time of updating dGCEWF_AvgE.
	double dGCEWF_AvgE, dGCEWF_SegAccumE;
	void ClusterDelt(unsigned int CTrial);
	void WolffClusterGCE(int Moves);//GCE cluster flipping with acceptance rates.
	//========== NVT Heat Path===========//
	int *QDeltaEi;//q=1,q=2...q=cQ
	//QVar
	bool SetHeatPath(double T, int *DEi, int n);
	int HeatPathDelt();
	void HeatPathTrial(int Moves);//NVT
	//========== N-Fold Way===========//

	//========== SE ===========//
	double * _Se_p;//store generalized S(e) function versus Energy
	bool SetSefun(double * Se_p, int size);
	void MetroplisTrialSe(int Moves);
	//========== WE ===========//
	double * _We_p;//store W(E)=exp[-S(E)] function versus Energy
	bool SetWefun(double * We_p, int size);
	void MetroplisTrialWe(int Moves);
	//======== Wang-Landau ==========//
	double * LnGe_p;//ln(g[e])
	unsigned long *GeN_p;//visited number of g[e]
	double Lnf, LnfPrecision, FlatRate;//percent
	//double dLowBound, dHighBound;

	void SetWL(double * LnGe_p, unsigned long *GeN_p, int size,
			double lnfPrecision, double flatRate);
	void WLResetGeN();
	bool WLCheckBelowAVG();
	void WangLandauTrial(int Moves);//Wand-Landau

	bool WLCheckBelowAVGBound();
	void WangLandauTrialBound(int Moves);
	void WarmupToBound();
	//======== Multicanonical ==========//
	double * Beta_p, *Gn_p;//_Se_p
	unsigned long *Hist_p;//visited number of g[e]
	double dh;//0<h<1; h -> 0
	unsigned int nFlatNonzero;
	int SetMul(double * beta_p, double *se_p, double *gn_p,
			unsigned long *hist_p, int size, double h);//len of beta_p,_Se_p,Hist_p,Gn_p =nEn,
	void MetroplisTrialSeBound(int Moves);//[nLowEBound, nHighEBound]
	void MulReSet();
	bool MulReweight();
};
#endif
