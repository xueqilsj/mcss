#include "mex.h"
#include "Ising.h"
#include "string.h"

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define ROUTE_IN prhs[0]
#define	ARG_IN1 prhs[1]
#define ARG_IN2 prhs[2]

#define ARG_OUT0 plhs[0]
#define ARG_OUT1 plhs[1]

static Ising *Ising_ptr = NULL;

#define RET_DOUBLE(_D_) \
    {\
	ARG_OUT0 = mxCreateDoubleMatrix(1, 1, mxREAL);\
    double *Arg = mxGetPr(ARG_OUT0);\
    Arg[0] = _D_;\
    }

void Cleanup(void) {
	/*>> clear Ising*/
#ifdef _DEBUG
	mexPrintf("Ising is terminating, destroying Ising_ptr\n");
#endif
	if (Ising_ptr != NULL)
		delete Ising_ptr;
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
		const mxArray *prhs[]) {

	//check route(get function name)
	if (mxCHAR_CLASS != mxGetClassID(ROUTE_IN))
		mexErrMsgTxt("Error: Invalid args");
	mwSize buflen = mxGetNumberOfElements(ROUTE_IN) + 1;
	char *buf = (char *) mxCalloc(buflen, sizeof(char));
	if (mxGetString(ROUTE_IN, buf, buflen) != 0) {
		mexErrMsgTxt("Could not convert string data.");
	}
	//require mwfree(buflen) at end.

	if (Ising_ptr == NULL) {
#ifdef _DEBUG
		mexPrintf("Ising is initializing, creating Ising_ptr\n");
#endif
		/* Create persistent Ising and register its cleanup. */
		Ising_ptr = new Ising();
		if (Ising_ptr == NULL)
			mexErrMsgTxt("Error: new Ising()");
		mexAtExit(Cleanup);
	}
#ifdef _DEBUG
	mexPrintf("Ising:N=%d*%d,hami=%g\n",Ising_ptr->nDim0,Ising_ptr->nDim1,
			Ising_ptr->Arg_io[Ising_ptr->AI_HAMI]);
#endif

	if (nrhs == 2) {
		if (strcmp(buf, "NVT") == 0) {
			int trials = mxGetScalar(ARG_IN1);
			Ising_ptr->MetroplisTrial(trials);
			RET_DOUBLE(Ising_ptr->Arg_io[Ising_ptr->AI_HAMI]);
		} else if (strcmp(buf, "SetTemperature") == 0) {
			Ising_ptr->SetTemperature(mxGetScalar(ARG_IN1));
		} else if (strcmp(buf, "InitConf") == 0) {
			Ising_ptr->InitConf(mxGetScalar(ARG_IN1));
		} else if (strcmp(buf, "new") == 0) {
			//Ising('new',conf)
			mwSize m, n;
			char *state_p = CheckConf(prhs, m, n);
			Ising_ptr->InitParameters();
			Ising_ptr->SetConf(state_p, m, n);
		} else {
			mexErrMsgTxt("input unknown");
		}
	} else if (nrhs == 1) {
		if (strcmp(buf, "Gauge") == 0) {
			Ising_ptr->Gauge();
		} else if (strcmp(buf, "GetArg_io") == 0) {
			ARG_OUT0 = mxCreateDoubleMatrix(1, Ising_ptr->AI_NUM, mxREAL);
			double *Arg = mxGetPr(ARG_OUT0);
			for (int i = 0; i < Ising_ptr->AI_NUM; i++)
				Arg[i] = Ising_ptr->Arg_io[i];
		} else if (strcmp(buf, "GetEnergyRange") == 0) {
			int dims[2];
			dims[0] = 1;
			dims[1] = Ising_ptr->nN + 1;
			ARG_OUT0 = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
			//TODO: int =32?
			int *Arg = (int *) mxGetPr(ARG_OUT0);
			for (int n = 0; n < dims[1]; n++) {
				Arg[n] = 4 * n - 2 * Ising_ptr->nN;//E(n)
			}

		} else {
			mexErrMsgTxt("input unknown");
		}
	} else {
		mexErrMsgTxt("input unknown");
	}
	mxFree(buf);//return normal
}
