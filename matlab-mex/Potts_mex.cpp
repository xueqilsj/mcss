#include "mex.h"
#include "Potts.h"
#include "string.h"

/*
 * for Matlab 2008
 sudo mv libstdc++.so.6 libstdc++.so.6.orig
 sudo ln -s /usr/lib/libstdc++.so.6.0.9 libstdc++.so.6
 sudo mv libgcc_s.so.1 libgcc_s.so.1.orig
 sudo ln -s /lib/libgcc_s.so.1

 >> mex -D_DEBUG potts.cpp Potts.cpp SeFun.cpp Random.cpp
 1000samples:0.3s new-delete every time

 */
#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define ROUTE_IN prhs[0]
#define	ARG_IN1 prhs[1]
#define ARG_IN2 prhs[2]

#define ARG_OUT0 plhs[0]
#define ARG_OUT1 plhs[1]

static Potts *Potts_ptr = NULL;

#define RET_DOUBLE(_D_) \
    {\
	ARG_OUT0 = mxCreateDoubleMatrix(1, 1, mxREAL);\
    double *Arg = mxGetPr(ARG_OUT0);\
    Arg[0] = _D_;\
    }

void Cleanup(void) {
	/*>> clear potts*/
#ifdef _DEBUG
	mexPrintf("Potts is terminating, destroying Potts_ptr\n");
#endif
	if (Potts_ptr != NULL)
		delete Potts_ptr;
}

char * CheckConf(const mxArray *prhs[], mwSize &m, mwSize &n) {
	if (mxGetElementSize(ARG_IN1) != 1) {
		mexErrMsgTxt("Element type of conf is int8");
	}
	if (mxIsClass(ARG_IN1, "sparse")) {
		mexErrMsgTxt("conf must be full matrix");
	}
	m = mxGetM(ARG_IN1);
	n = mxGetN(ARG_IN1);
	if (MIN(m,n) < 3) {
		mexErrMsgTxt("every dimension of conf be greater than 3");
	}
	return (char *) mxGetPr(ARG_IN1);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	//check route(get function name)
	if (mxCHAR_CLASS != mxGetClassID(ROUTE_IN))
		mexErrMsgTxt("Error: Invalid args");
	mwSize buflen = mxGetNumberOfElements(ROUTE_IN) + 1;
	char *buf = (char *) mxCalloc(buflen, sizeof(char));
	if (mxGetString(ROUTE_IN, buf, buflen) != 0) {
		mexErrMsgTxt("Could not convert string data.");
	}
	//require mwfree(buflen) at end.

	if (Potts_ptr == NULL) {
#ifdef _DEBUG
		mexPrintf("potts is initializing, creating Potts_ptr\n");
#endif
		/* Create persistent Potts and register its cleanup. */
		Potts_ptr = new Potts();
		if (Potts_ptr == NULL)
			mexErrMsgTxt("Error: new Potts()");
		mexAtExit(Cleanup);
	}
#ifdef _DEBUG
	mexPrintf("Potts:N=%d*%d,Q=%d,hami=%g\n",Potts_ptr->nDim0,Potts_ptr->nDim1,
			int(Potts_ptr->cQ), Potts_ptr->Arg_io[Potts_ptr->AI_HAMI]);
#endif

	if (nrhs == 2) {
		if (strcmp(buf, "FAFH") == 0) {
			int trials = mxGetScalar(ARG_IN1);
			Potts_ptr->MetroplisTrialSe(trials);
			RET_DOUBLE(Potts_ptr->Arg_io[Potts_ptr->AI_HAMI]);
		} else if (strcmp(buf, "WL") == 0) {
			int trials = mxGetScalar(ARG_IN1);
			Potts_ptr->WangLandauTrial(trials);
			RET_DOUBLE(Potts_ptr->Arg_io[Potts_ptr->AI_HAMI]);
		} else if (strcmp(buf, "NVT") == 0) {
			int trials = mxGetScalar(ARG_IN1);
			Potts_ptr->MetroplisTrial(trials);
			RET_DOUBLE(Potts_ptr->Arg_io[Potts_ptr->AI_HAMI]);
		} else if (strcmp(buf, "SetTemperature") == 0) {
			Potts_ptr->SetTemperature(mxGetScalar(ARG_IN1));
		} else if (strcmp(buf, "InitConf") == 0) {
			Potts_ptr->InitConf(mxGetScalar(ARG_IN1));
		} else if (strcmp(buf, "WLUpdatePrecision") == 0) {
			RET_DOUBLE(Potts_ptr->Lnf);
			Potts_ptr->Lnf = mxGetScalar(ARG_IN1);
		} else {
			mexErrMsgTxt("input unknown");
		}
	} else if (nrhs == 1) {
		if (strcmp(buf, "WLCheckBelowAVG") == 0) {
			int dims[2] = { 1, 1 };
			ARG_OUT0 = mxCreateNumericArray(2, dims, mxLOGICAL_CLASS, mxREAL);
			bool *below = (bool*) mxGetData(ARG_OUT0);
			below[0] = Potts_ptr->WLCheckBelowAVG();
		} else if (strcmp(buf, "GetArg_io") == 0) {
			ARG_OUT0 = mxCreateDoubleMatrix(1, Potts_ptr->AI_NUM, mxREAL);
			double *Arg = mxGetPr(ARG_OUT0);
			for (int i = 0; i < Potts_ptr->AI_NUM; i++)
				Arg[i] = Potts_ptr->Arg_io[i];
		} else if (strcmp(buf, "GetEnergyRange") == 0) {
			int dims[2];
			dims[0] = 1;
			dims[1] = Potts_ptr->nN + Potts_ptr->nN + 1;
			ARG_OUT0 = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
			//TODO: int =32?
			int *Arg = (int *) mxGetPr(ARG_OUT0);
			for (int n = 0; n < dims[1]; n++) {
				Arg[n] = n - dims[1] + 1;//E(n)
			}
		} else if (strcmp(buf, "Gauge") == 0) {
			Potts_ptr->Gauge();
		} else if (strcmp(buf, "WLResetGeN") == 0) {
			Potts_ptr->WLResetGeN();
		} else {
			mexErrMsgTxt("input unknown");
		}
	} else if (nrhs = 3) {
		if (strcmp(buf, "new") == 0) {
			//potts('new',conf,Q)
			mwSize m, n;
			char *state_p = CheckConf(prhs, m, n);
			char Q = mxGetScalar(ARG_IN2);
			Potts_ptr->InitParameters(Q);
			Potts_ptr->SetConf(state_p, m, n);
		} else if (strcmp(buf, "SetWl") == 0) {
			//('SetWL',lnfPrecision, flatRate)
			if (nlhs != 2)
				mexErrMsgTxt("lnge,gen output");
			int dims[2];
			dims[0] = 1;
			dims[1] = Potts_ptr->nN + Potts_ptr->nN + 1;
			ARG_OUT0 = mxCreateDoubleMatrix(1, dims[1], mxREAL);
			ARG_OUT1 = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
			Potts_ptr->SetWL((double *) mxGetPr(ARG_OUT0),
					(unsigned long *) mxGetPr(ARG_OUT1), dims[1], mxGetScalar(
							prhs[1]), mxGetScalar(prhs[2]));
#ifdef _DEBUG
			mexPrintf("size,ulong=%d,uint=%d\n", sizeof(unsigned long),
					sizeof(unsigned int));
#endif
		} else {
			mexErrMsgTxt("input unknown");
		}
	} else {
		mexErrMsgTxt("input unknown");
	}
	mxFree(buf);//return normal
}
