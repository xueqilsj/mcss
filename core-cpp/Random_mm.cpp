#include  <cstdlib>
#include  <ctime>
#include  <climits> //INT_MAX
#include  <cmath>
#include  "Random_mm.h"

Random::Random() {
	Srand((unsigned) time(NULL));
}
/// set seed
void Random::Srand(uint32_t seed) {
	int n;
	for (n = 0; n < 55; n++)
		ia[n] = conv1 * (seed = (a * seed + c) % m);
	p = 0;
	pp = 24;
}
double Random::Real() {
	if (--p < 0)
		p = 54;
	if (--pp < 0)
		pp = 54;
	return conv2 * (ia[p] += ia[pp]);
}
/// return a double random number
double Random::Real(double dfLow, double dfHigh) {
	if (dfLow > dfHigh)
		Real(dfHigh, dfLow);

	double ran = INT_MAX;///(1<<(sizeof(int)*8-1))-1;[0,(2^31 - 1)]
	ran = ((dfHigh - dfLow) * rand()) / ran + dfLow;

	return ran;
}
/// return a integer random number
int32_t Random::Number(int32_t nLow, int32_t nHigh) {
	if (--p < 0)
		p = 54;
	if (--pp < 0)
		pp = 54;
	return (ia[p] += ia[pp]) % (nHigh - nLow + 1) + nLow;
}

