/*
 * RandWalker.cpp
 *  Basic C++ Class for Random walker model
 *  Created on: May 20, 2009
 *
 * remain precision:
 * double a=atof("-1.42343454566");
 * double b;
 * sscanf("-1.42343454566","%lf",&b);
 *
 */
#include "RandWalker.h"
#include "nvector.h"

#include <iostream> //< for cout
#include <cmath>
#include <ctype.h> //isdigit
using namespace std;
RandWalker::RandWalker() {

}
void RandWalker::SetConf(char * State_p, int dim0, int dim1) {
	nDim0 = dim0;
	nDim1 = dim1;
	nN = nDim0 * nDim1;
	States_p = State_p;//delete States_p before assign new one.
}
void RandWalker::InitConf(unsigned int *Snake_pos, unsigned int len, int d) {//d<0 lower state;d==0 custom;d>0 upper state
	if (States_p == NULL) {
		cout << "#Error: Set configures first by \'SetConf()\'" << endl;
		return;
	}
	if (Snake_pos == NULL || len > nDim1) {
		cout << "#Error: InitConf input args" << endl;
		return;
	}

	int dim = nDim0 / 2;
	if (d < 0) {
		for (unsigned int l = 0; l < len; l++)
			Snake_pos[l] = l;
	} else if (d > 0) {
		for (unsigned int l = 0; l < len; l++)
			Snake_pos[l] = dim * nDim1 + l;
	}
	for (unsigned int n = 0; n < nN; n++) {
		States_p[n] = BackGround_V;
	}
	SnakePos_p = Snake_pos;
	nInitconf = d;
	nLen = len;
	nSnakeHead_i = nLen - 1;

	for (unsigned int n = 0; n < len; n++) {
		States_p[SnakePos_p[n]] = Face_V;
	}

}
int RandWalker::MetroplisTrial(int Moves) {
	int CurStep = 0;
	unsigned int acn = 0;
	unsigned int neib;
	int go_i;
	while (CurStep < Moves) {

		if (nSnakeHead_i == nLen)
			nSnakeHead_i -= nLen;

		acn = 0;
		neib = SnakePos_p[nSnakeHead_i] + 1;//right
		if (neib % nDim1 == 0)
			neib -= nDim1;

		if (States_p[neib] == BackGround_V)
			Acc_dir[acn++] = neib;

		neib = SnakePos_p[nSnakeHead_i] - 1;//left
		if (nHead % nDim1 == 0)
			neib += nDim1;
		if (States_p[neib] == BackGround_V)
			Acc_dir[acn++] = neib;

		neib = SnakePos_p[nSnakeHead_i] - nDim1;//up
		if (nHead < nDim1)
			neib += nN;
		if (States_p[neib] == BackGround_V)
			Acc_dir[acn++] = neib;

		neib = SnakePos_p[nSnakeHead_i] + nDim1;//down
		if (neib >= nN)
			neib -= nN;
		if (States_p[neib] == BackGround_V)
			Acc_dir[acn++] = neib;

		if (acn == 0)
		{
			return acn;
		}
		else (acn==1)
		{
			go_i=0;
		}
		else
		{
			go_i = ran.Number(0, acn - 1);
		}

		States_p[Acc_dir[go_i]] = Face_V;
		if (nSnakeHead_i == 0) {
			States_p[SnakePos_p[nLen - 1]] = BackGround_V;
			SnakePos_p[nLen - 1] = Acc_dir[go_i];
		} else {
			States_p[SnakePos_p[nSnakeHead_i - 1]] = BackGround_V;
			SnakePos_p[nSnakeHead_i - 1] = Acc_dir[go_i];
		}

		nSnakeHead_i++;

		uCurMove++;
		CurStep += 1;
	}
	return acn;
}//index exp(-dE)
RandWalker::~RandWalker() {
}

