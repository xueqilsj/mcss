#include  <cstdlib>
#include  <ctime>
#include  <climits> //INT_MAX
#include  <cmath>
#include  "Random_cl.h"

Random::Random() {

	srand((unsigned) time(NULL));

}
/// set seed
void Random::Srand(unsigned seed) {
	if (seed > 1) {
		srand(seed);
	} else if (seed == 1) {
		srand((unsigned) time(NULL));
	}
}
double Random::Real() {
	return rand() / double(INT_MAX);
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
int Random::Number(int nLow, int nHigh) {
	if (nLow > nHigh)
		Number(nHigh, nLow);

	return rand() % (nHigh - nLow + 1) + nLow;
}


