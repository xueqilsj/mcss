/*
 ============================================================================
 Name        : NVTSPart.cpp
 Description : Sampling in NVT Single Particle model with Metroplis Algorithm
 ============================================================================
 */

#include "SParticle.h"
#include "nvector.h"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

using namespace std;

class CNVTSPart: public SParticle {

public:
	char filename[128];
	CNVTSPart(double initPt, double off, double T) {
		InitParameters(initPt, off, T);
	}
	~ CNVTSPart() {

	}
	void RunMetroplis_Pak(unsigned long Samples, unsigned int Stride) {
		Gauge_Pak();//initp= 0, offset=0.002, temp=1/5.0
		sprintf(filename, "nvtspak%lg_%lg~%lg_%ld_%dTraj.txt", dX, dOffset,
				Arg_io[AI_TEMPER], Samples, Stride);

		std::ofstream ofile(filename);
		for (unsigned long CurSample = 0; CurSample < Samples; CurSample++) {
			MetroplisTrial_Pak(Stride);
			ofile << setprecision(16) << Arg_io[AI_HAMI] << '\t' << dX << endl;
		}
		ofile.close();

	}
	void RunMetroplis_MP(unsigned long Samples, unsigned int Stride) {
		Gauge_MP();//initp= -1.3, offset=0.1, temp=0.4
		sprintf(filename, "nvtsmp%lg_%lg~%lg_%ld_%dTraj.txt", dX, dOffset,
				Arg_io[AI_TEMPER], Samples, Stride);

		std::ofstream ofile(filename);
		for (unsigned long CurSample = 0; CurSample < Samples; CurSample++) {
			MetroplisTrial_MP(Stride);
			ofile << setprecision(16) << Arg_io[AI_HAMI] << '\t' << dX << endl;
		}
		ofile.close();

	}

};
int main(int argc, char *argv[]) {
	if (argc == 7) {
		clock_t __start = clock();
		CNVTSPart p(atof(argv[1]), atof(argv[2]), atof(argv[3]));

		if (argv[6][0] == 'p') {
			p.RunMetroplis_Pak(atoi(argv[4]), atoi(argv[5]));
		} else if (argv[6][0] == 'm') {
			p.RunMetroplis_MP(atoi(argv[4]), atoi(argv[5]));
		} else {
			cout << "unknow flag: " << argv[5] << endl;
		}
		cout << ((double) (clock() - __start)) / CLOCKS_PER_SEC << "\t#"
				<< p.filename << std::endl;

	} else {
		cout << "<NVTSPart>  initX offset T samples stride m/p" << endl;//stride =MCSweep/site
		return -1;
	}
	return 0;
}

