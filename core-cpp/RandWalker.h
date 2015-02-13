/*
 * RandWalker.h
 *
 *  Created on: Dec 20, 2009
 *
 *   0---N-->
 * (0,0)----y,dim1--->
 *   |
 *   |
 *   x,
 *  dim0
 *   |
 *   v
 */

#ifndef _RandWalker_H
#define _RandWalker_H

#ifdef RANDOM_MM
#include "Random_mm.h"
#elif defined(RANDOM_CL)
#include "Random_cl.h"
#else
#include "Random.h"
#endif

#define ACC_MAXDIR 3
class RandWalker {
protected:

public:

	double Acc_dir[ACC_MAXDIR + 1];
	enum {
		BackGround_V = 1, Face_V = 2
	};
	unsigned int nDim0, nDim1, nN, nCurTrial;// Dimensions and size
	char * States_p;//configuration, size=nDim0*nDim1;
	char * SnakePos_p;//configuration, size=nDim0*nDim1;

	unsigned int nSnakeHead_i, nLen;

	int nInitconf;//starting configuration:
	unsigned long long uCurMove;
	int nCurHami, nDeltaE;//nEn:the number of total energy levels

	Random ran;//random number generator,automatic initiation.

	RandWalker();
	virtual ~RandWalker();

	//void InitParameters(char q);

	void SetConf(char * State_p, int dim0, int dim1);
	void InitConf(unsigned int *Snake_pos, int len, int d);//d<0 lower state;d==0 custom;d>0 upper state


	void MetroplisTrial(int Moves);//index exp(-dE)

};
#endif
