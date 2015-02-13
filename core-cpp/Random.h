#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <limits.h>

#ifdef _WIN32
typedef unsigned int uint32_t;
typedef signed int int32_t;
#else
#include <stdint.h> // for uint32_t
#endif

class Random {
private:
	/* Period parameters */
	static const uint32_t N = 624;
	static const int32_t M = 397;
	static const uint32_t MATRIX_A = (uint32_t) 0x9908b0df; /* constant vector a */
	static const uint32_t UPPER_MASK = (uint32_t) 0x80000000; /* most significant w-r bits */
	static const uint32_t LOWER_MASK = (uint32_t) 0x7fffffff; /* least significant r bits */
	uint32_t mag01[2];
	uint32_t mt[N]; /* the array for the state vector  */
	uint32_t mti; /* mti==N+1 meansmt[N] is not initialized*/
public:
	uint32_t Seed;
	Random(); /* initializes mt[N] with a seed */
	void Srand(uint32_t s);
	// Int31: a random number on [0,0x7fffffff]-Interval
	int32_t Number(int32_t nLow, int32_t nHigh);
	//a random number on [0,1]-real-Interval */
	double Real();
	double Real(double dfLow, double dfHigh);//poor!!
	//"Unsigned Int 32bit on [0,0xffffffff]-Interval */
	uint32_t uNumber32();
};
#endif
