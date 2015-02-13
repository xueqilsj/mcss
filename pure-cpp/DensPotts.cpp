/*
 ============================================================================
 Name        : GCEPotts.cpp
 Description : To enumerate the configurations of 2D Q-state Potts model.

 L*L=4*4,Q=4,	 237.77sec.	#densp4_4_4His.txt
 L*L=5*5,Q=4,	 20days 	#densp5_4_4His.txt
 ============================================================================
 */

#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "Potts.h"
#include "nvector.h"

using namespace std;

class DensPotts: public Potts {
protected:
	unsigned long long *HisE;
public:
	char filename[128];
	DensPotts() {
	}
	DensPotts(int dim0, int dim1, char Q) {
		InitParameters(Q);
		States_p = new char[dim0 * dim1 + 1];
		SetConf(States_p, dim0, dim1);
		InitConf(-1);//set all spin to q=1
	}
	~ DensPotts() {
		delete[] States_p;
	}
	void CalSpins(int curE, unsigned int start) {
		//conf: 1/2 1 1 1 1 ... 1 1
		if (start == 0) {
			HisE[nEn + curE - 1]++;
			for (char q = 2; q <= cQ; q++) {
				QVar = q;
				DeltS(start);
				States_p[start] = q;
				curE += nDeltaE;
				HisE[nEn + curE - 1]++;
			}
			States_p[start] = 1;//recover

		} else {
			CalSpins(curE, start - 1);
			for (char q = 2; q <= cQ; q++) {
				QVar = q;
				DeltS(start);
				States_p[start] = q;
				curE += nDeltaE;
				CalSpins(curE, start - 1);
				//Gauge();
				//cout << '<' << curE << ':' << nCurHami << '>' << endl;
			}
			States_p[start] = 1;//recover
		}
	}
	void RunCalDens() {
		time_t first = time(NULL);
		HisE = new unsigned long long[nEn];
		for (int n = 0; n < nEn; n++) {
			HisE[n] = 0;
		}
		/*for (int n = 0; n < nN; n++) {
		 States_p[n] = 1;
		 }//InitConf(-1);//set all spin to q=1
		 */
		Gauge();
		CalSpins(nCurHami, nN - 1);
		nvector<unsigned long long> his(HisE, nEn);

		double sec = difftime(time(NULL), first);
		sprintf(filename, "densp%d_%d_%d_%lgHis.txt", nDim0, nDim1, int(cQ),
				sec);
		his.dump(filename);
		cout << his << endl;
		delete[] HisE;

	}

};

int main(int argc, char *argv[]) {

	int d0, d1, q;
	if (argc == 4) {
		d0 = atoi(argv[1]);
		d1 = atoi(argv[2]);
		q = atoi(argv[3]);
		if (d0 < 3 || d1 < 3) {
			cout << "Dimensions too small: " << d0 << '*' << d1 << endl;
			return -1;
		} else if (q < 2) {
			cout << "Q<2: " << q << endl;
			return -1;
		}
		DensPotts p(d0, d1, q);
		p.RunCalDens();
	} else {
		cout << "<DensPotts> Dim0 Dim1 Q" << endl;//sampling at first time
		return -1;
	}
	return 0;
}

