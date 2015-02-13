#ifndef _RANDOM_MM_H_
#define _RANDOM_MM_H_

// additive lagged Fibonacci generator of mitchell and Moore,described by knuth.

#include <stdint.h> // for uint32_t
class Random {
private:
	/* Period parameters */
	static const uint32_t a = 2416;
	static const uint32_t c = 374441;
	static const uint32_t m = 1771875;
	static const double conv1 = 2423.9674;
	static const double conv2 = 1 / 4294967296.0;
	uint32_t ia[55];
	long p, pp;
public:
	Random();
	void Srand(uint32_t seed); //1: let seed=time(NULL); >1: using the given seed
	double Real(double dfLow, double dfHigh);//< inclusive boundaries
	double Real();//< inclusive boundaries
	int32_t Number(int32_t nLow = 0, int32_t nHigh = 255);//< inclusive boundaries
};

#endif
