/*
 * SeFun.h
 *
 *  Created on: Mar 9, 2009
 */

#ifndef SeFun_H_
#define SeFun_H_

#include <float.h>
#include <vector>
#include <iterator>

#define DOUBLE_EQ(x,v) (((v - DBL_EPSILON) < x) && (x <( v + DBL_EPSILON)))
#define DELTA_ENERGY 1.0

struct TEP {
public:
	double T;
	double E;
	double Sigma;
	double A;
	double SE;
};

class SeFun {
protected:
	std::vector<TEP> _TELine_l;//low->high
public:
	int _nExDir;
	int _uSizeN;
	double dLflatEx, dHflatEx;
	bool bSmooth;
	SeFun();
	virtual ~SeFun();
	void ReSetFun(int nNum0, int nNum1, double *dTE);
	void InitFun(double lflatEx, double hflatEx, bool smooth);
	double SEFunction(double energy, double *beta, double High_t = 0.89,
			double Low_t = 0.62);
	//double SEFunction(double curEnergy, double deltaEnergy);
	double _DIntegralBeta(std::vector<TEP>::iterator begin,
			std::vector<TEP>::iterator end, double e, double *beta);
};
#endif /* SeFun_H_ */

