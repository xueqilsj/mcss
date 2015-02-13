#include "mex.h"
#include "Linear.h"
#include "string.h"

#define ROUTE_IN prhs[0]
#define	ARG_IN1 prhs[1]
#define ARG_IN2 prhs[2]

#define ARG_OUT0 plhs[0]
#define ARG_OUT1 plhs[1]

static Linear *Linear_ptr = NULL;

#define RET_DOUBLE(_D_) \
    {\
	ARG_OUT0 = mxCreateDoubleMatrix(1, 1, mxREAL);\
    double *Arg = mxGetPr(ARG_OUT0);\
    Arg[0] = _D_;\
    }

void Cleanup(void) {
	/*>> clear Linear*/
#ifdef _DEBUG
	mexPrintf("Linear is terminating, destroying Linear_ptr\n");
#endif
	if (Linear_ptr != NULL)
		delete Linear_ptr;
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

	if (Linear_ptr == NULL) {
#ifdef _DEBUG
		mexPrintf("Linear is initializing, creating Linear_ptr\n");
#endif
		/* Create persistent Linear and register its cleanup. */
		Linear_ptr = new Linear();
		if (Linear_ptr == NULL)
			mexErrMsgTxt("Error: new Linear()");
		mexAtExit(Cleanup);
	}
#ifdef _DEBUG
	mexPrintf("Linear:dim=%g~%g,offset=%g,hami=%g\n",Linear_ptr->dDim[0],Linear_ptr->dDim[1],
			Linear_ptr->dOffset[0], Linear_ptr->Arg_io[Linear_ptr->AI_HAMI]);
#endif

	if (nrhs == 2) {
		//Linear('type',trials)
		if (strcmp(buf, "NVT") == 0) {
			int trials = mxGetScalar(ARG_IN1);
			Linear_ptr->MetroplisTrial(trials);
			RET_DOUBLE(Linear_ptr->Arg_io[Linear_ptr->AI_HAMI]);
		} else if (strcmp(buf, "GetXY") == 0) {
			if (nlhs != 2)
				mexErrMsgTxt("output [X,Y]");
			double step = mxGetScalar(ARG_IN1);
			int dims[2];
			dims[0] = 1;
			dims[1] = (Linear_ptr->dDim[1] - Linear_ptr->dDim[0]) / step;
			ARG_OUT0 = mxCreateDoubleMatrix(1, dims[1], mxREAL);
			ARG_OUT1 = mxCreateDoubleMatrix(1, dims[1], mxREAL);
			double *x = (double *) mxGetPr(ARG_OUT0);
			double *y = (double *) mxGetPr(ARG_OUT1);
			double old = Linear_ptr->dPoint[0];
			for (int n = 0; n < dims[1]; n++) {
				Linear_ptr->dPoint[0] = Linear_ptr->dDim[0] + n * step;
				Linear_ptr->Gauge();
				y[n] = Linear_ptr->Arg_io[Linear_ptr->AI_HAMI];
				x[n] = Linear_ptr->dPoint[0];
			}
			Linear_ptr->dPoint[0] = old;
			Linear_ptr->Gauge();
		} else if (strcmp(buf, "GetPoint") == 0) {
			RET_DOUBLE(Linear_ptr->GetPoint(mxGetScalar(ARG_IN1)));
		}
	} else if (nrhs == 1) {
		if (strcmp(buf, "Gauge") == 0) {
			Linear_ptr->Gauge();
		} else if (strcmp(buf, "Delt") == 0) {
			Linear_ptr->Delt();
		} else if (strcmp(buf, "Accepted") == 0) {
			Linear_ptr->Accepted();
		} else if (strcmp(buf, "GetArg_io") == 0) {
			ARG_OUT0 = mxCreateDoubleMatrix(1, Linear_ptr->AI_NUM, mxREAL);
			double *Arg = mxGetPr(ARG_OUT0);
			for (int i = 0; i < Linear_ptr->AI_NUM; i++)
				Arg[i] = Linear_ptr->Arg_io[i];
		} else {
			mexErrMsgTxt("input unknown");//mxFree(buf) by Matlab;
		}
	} else if (nrhs = 4) {
		if (strcmp(buf, "InitParameters") == 0) {
			//Linear('InitParameters', initp, offset, temp)
			Linear_ptr->InitParameters(mxGetScalar(prhs[1]), mxGetScalar(prhs[2]),
					mxGetScalar(prhs[3]));
		} else {
			mexErrMsgTxt("input unknown");
		}
	} else {
		mexErrMsgTxt("input unknown");
	}
	mxFree(buf);//return normal
}
