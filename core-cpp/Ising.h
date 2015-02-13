/*
 * Ising.h
 *
 *  Created on: Jun 2, 2009
 *      Hamitonian for Ising model
 *      H=-j \sigma_<ij> s_i*s_j-B\sigma_i s_i
 */

#ifndef _Ising_H
#define _Ising_H

#ifdef RANDOM_MM
#include "Random_mm.h"
#elif defined(RANDOM_CL)
#include "Random_cl.h"
#else
#include "Random.h"
#endif

#define TO32       4294967296.0  //2^32
//J=1
class Ising {
protected:
	double _dDeltaExp[5];//-8,-4,0,4,8
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
	unsigned int nDim0, nDim1, nN, nCurTrial;
	int nEn, nCurHami, nCurMag, nAccepted;//< Dimension and size
	//nN:size of model,nEn: the total number of energy levels
	double dBeta, dAlpha, dUn;//for Generalized canonical ensemble (dUn=U/N).
	unsigned long long uCurMove;
	int nLowEBound, nHighEBound;//[nLowEBound, nHighEBound],for flat histogram:WL,Mul
	int nDeltaE, nDeltaMag;
	char * States_p;//< external point
	int nInitconf;
	double Arg_io[AI_NUM];
	Random ran, ran2;//< random number generator
	Ising();
	virtual ~Ising();

	void SetArg(int site, double df);
	double GetArg(int site);

	void InitParameters();
	void SetConf(char * State_p, int dim0, int dim1);
	void InitConf(int d);//d<0 ground state;d==0 random state;d>0 anti-ground state
	bool FromConf(const char * ConfName);
	bool DumpConf(const char * ConfName);
	bool LoadConf(const char * ConfName);
	void Gauge();
	int Delt();
	void DeltS(unsigned int CTrial);//for Enumerating Conf.
	//========== NVT ===========//
	void SetTemperature(double T);
	void MetroplisTrial(int Moves);//NVT
	//========== GCE Metroplis===========//
	void MetroplisTrialGCE(int Moves);
	//========== NVT Heat Bath===========//
	unsigned int uLocalFieldPro[5];//[-4,-2,0,2,4]
	void SetHeatBath(double T);
	void HeatBathTrial(int Moves);//NVT
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

	void WLResetGeNBound();
	bool WLCheckBelowAVGBound();
	void WangLandauTrialBound(int Moves);
	void WarmupToBound();
	//===========Wolff Cluster=======//
	int * ClusterStack_p;//<
	double dPadd;
	//	int downBonds, upBonds;
	void WolffCluster(int Moves);//NVT

	int * ClusterSeq_p;//roll back if GCE cluster flipping refused.

	int csp;//the size of current cluster
	int nSegAccumCS;//accumulate size a part of cluster
	double dAccumCS;

	int nGCEWF_AvgC;//the time of updating dGCEWF_AvgE.
	double dGCEWF_AvgE, dGCEWF_SegAccumE;
	void ClusterDelt(unsigned int CTrial);
	void WolffClusterGCE(int Moves);//GCE cluster flipping with acceptance rates.
};

#endif
