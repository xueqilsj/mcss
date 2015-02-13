/*
 A program for MT19937, with initialization improved 2002/1/26.
 Coded by Takuji Nishimura and Makoto Matsumoto.
 */
#include "Random.h"
#include <ctime>
#ifndef INT32_T_MAX
#define INT32_T_MAX 2147483647
#endif

Random::Random() {
	mti = N + 1;
	Seed = time(NULL);
	Srand(Seed);

	mag01[0] = 0;
	mag01[1] = MATRIX_A;
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

}
/* initializes mt[N] with a seed */
void Random::Srand(uint32_t s) {
	Seed = s;
	mt[0] = Seed;

	for (mti = 1; mti < N; mti++) {
		mt[mti] = ((uint32_t) 1812433253 * (mt[mti - 1] ^ (mt[mti - 1] >> 30))
				+ mti);
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].       */
		/* 2002/01/09 modified by Makoto Matsumoto             */
	}
}

// Int31: a random number on [0,0x7fffffff]-Interval
int32_t Random::Number(int32_t nLow, int32_t nHigh) {
	//if (nLow > nHigh) {Number(nHigh, nLow);}
	//int32_t ran = uNumber32() >> 1;
	return (uNumber32() >> 1) % (nHigh - nLow + 1) + nLow;
}

//a random number on [0,1]-real-Interval */
double Random::Real() {
	return uNumber32() * ((double) 1.0 / 4294967295.0);
	/* divided by 2^32-1 */
}

double Random::Real(double dfLow, double dfHigh)//< inclusive boundaries
{
	return ((dfHigh - dfLow) * Real()) + dfLow;
}

//"Unsigned Int 32bit on [0,0xffffffff]-Interval */
uint32_t Random::uNumber32() {
	uint32_t y;

	if (mti >= N) { /* generate N words at one time */
		uint32_t kk;

		if (mti == N + 1) /* if init_genrand() has not been called, */
			Srand(5489); /* a default initial seed is used */

		for (kk = 0; kk < N - M; kk++) {
			y = ((mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK)) >> 1;
			mt[kk] = mt[kk + M] ^ mag01[mt[kk + 1] & 1] ^ y;
		}
		for (; kk < N - 1; kk++) {
			y = ((mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK)) >> 1;
			mt[kk] = mt[kk + (M - N)] ^ mag01[mt[kk + 1] & 1] ^ y;
		}
		y = ((mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK)) >> 1;
		mt[N - 1] = mt[M - 1] ^ mag01[mt[0] & 1] ^ y;

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680;
	y ^= (y << 15) & 0xefc60000;
	y ^= (y >> 18);

	return y;
}

