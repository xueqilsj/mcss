/*
 ============================================================================
 Name        : Densising.cpp
 Description : To enumerate the configurations of 2D Ising model.

 ============================================================================
 */

#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "Ising.h"
#include "nvector.h"

using namespace std;

class DensIsing: public Ising {
protected:
	unsigned long long *HisE;
public:
	char filename[128];
	DensIsing() {
	}
	DensIsing(int dim0, int dim1) {
		InitParameters();
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(-1);//set all spin to q=1
	}
	~ DensIsing() {
		delete[] States_p;
	}
	void CalSpins(int curE, unsigned int start) {
		if (start == 0) {
			HisE[(2 * nN + curE) / 4]++;
			DeltS(start);
			curE += nDeltaE;
			HisE[(2 * nN + curE) / 4]++;
		} else {
			CalSpins(curE, start - 1);
			DeltS(start);
			curE += nDeltaE;
			States_p[start] *= -1;
			CalSpins(curE, start - 1);
			//Gauge();
			//cout << '<' << curE << ':' << nCurHami << '>' << endl;
			States_p[start] *= -1;
		}
	}
	void RunCalDens() {
		time_t first = time(NULL);
		HisE = new unsigned long long[nEn];
		for (int n = 0; n < nEn; n++) {
			HisE[n] = 0;
		}
		Gauge();
		CalSpins(nCurHami, nN - 1);
		nvector<unsigned long long> his(HisE, nEn);

		double sec = difftime(time(NULL), first);
		sprintf(filename, "densi%d_%d_%lgHis.txt", nDim0, nDim1, sec);
		his.dump(filename);
		cout << his << endl;
		delete[] HisE;

	}

};

int main(int argc, char *argv[]) {
	int d0, d1;
	if (argc == 3) {
		d0 = atoi(argv[1]);
		d1 = atoi(argv[2]);
		if (d0 < 3 || d1 < 3) {
			cout << "Dimensions too small: " << d0 << '*' << d1 << endl;
			return -1;
		}
		DensIsing p(d0, d1);
		p.RunCalDens();
	} else {
		cout << "<DensIsing> Dim0 Dim1" << endl;//sampling at first time
		return -1;
	}
	return 0;
}

