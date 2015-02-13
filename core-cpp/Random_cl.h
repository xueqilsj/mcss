#ifndef _RANDOM_CL_H_
#define _RANDOM_CL_H_

//Random based on C library

class Random {

public:
	Random();
	void Srand(unsigned seed); //1: let seed=time(NULL); >1: using the given seed
	double Real(double dfLow, double dfHigh);//< inclusive boundaries
	double Real();//< inclusive boundaries
	int Number(int nLow = 0, int nHigh = 255);//< inclusive boundaries
};

#endif
