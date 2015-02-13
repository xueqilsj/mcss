%module (docstring="This is a Python SWIG-wrapped Monte Carlo module, Jun 20,2009") mcss
%{
#include "Random.h"
#include <Python.h>
#include <numpy/arrayobject.h>   /* numpy  as seen from C */
#include "Potts.h"
#include "Ising.h"
#include "SParticle.h"
#include "Statis.h"
#include "LqFluid.h"
#include <limits>
#include <iostream>
%}

%apply double *OUTPUT { double *DeltHami,double *beta };//before %include
%include "Random.h"
//%include "SeFun.h"
%include "Potts.h"
%include "Ising.h"
%include "SParticle.h"
%include "Statis.h"
%include "LqFluid.h"

#define QUOTE(s) # s      /* turn s into string "s" */
#define NDIM_CHECK(a, expected_ndim) \
       if (a->nd != expected_ndim) { \
          PyErr_Format(PyExc_ValueError, \
          "%s array is %d-dimensional, expected to be %d-dimensional",\
                         QUOTE(a), a->nd, expected_ndim); \
          return NULL; \
       }
#define DIM_CHECK(a, dim, expected_length) \
       if (a->dimensions[dim] != expected_length) { \
          PyErr_Format(PyExc_ValueError, \
          "%s array has wrong %d-dimension=%d (expected %d)", \
                     QUOTE(a),dim,a->dimensions[dim],expected_length); \
          return NULL; \
       }
#define DIM_GT_CHECK(a, dim, expected_length) \
       if (a->dimensions[dim] < expected_length) { \
          PyErr_Format(PyExc_ValueError, \
          "%s array has wrong %d-dimension=%d (expected %d)", \
                     QUOTE(a),dim,a->dimensions[dim],expected_length); \
          return NULL; \
       }
#define TYPE_CHECK(a, tp) \
       if (a->descr->type_num != tp) { \
          PyErr_Format(PyExc_TypeError, \
          "%s array is not of correct type (%d)", QUOTE(a), tp); \
          return NULL; \
       }
#define NULL_CHECK(a,message) \
		if (a==NULL) {\
			PyErr_Format(PyExc_TypeError, \
        	"%s is NULL, %s", QUOTE(a),message); \
          	return NULL;\
     	}
     	
%extend Potts
{
	PyObject* PySetConf(PyObject* a_)
	{
		PyArrayObject* a = (PyArrayObject*) PyArray_ContiguousFromObject(a_, NPY_INT8, 0, 0);
		NULL_CHECK(a,"PySetConf:PyArray_ContiguousFromObject(numpy.int8)");
		NDIM_CHECK(a,2);
		self->SetConf(a->data,a->dimensions[0],a->dimensions[1]);
		return PyArray_Return(a);
	}
	PyObject* PyGetQDis()
	{
		PyArrayObject *QDis; 
		npy_intp dims[1];
		dims[0]  = self->cQ + 1;
		QDis = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_UINT);
		unsigned int * x=(unsigned int *)QDis->data;
		double QM2 = self->QDistri(x, self->cQ + 1);
		double QM1 = self->cQ * x[0];
		QM1=(QM1 / self->nN - 1) / (self->cQ - 1);
		return Py_BuildValue("(ddO)",QM1,QM2,QDis); 
	}
				
	PyObject* PyGetEnergyRange()
	{
		PyArrayObject *e; npy_intp dims[1];
		dims[0]  = self->nN+self->nN+1;
		e = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
		int *x;
		for(int n=0;n<dims[0];n++)
		{
			x=(int *)e->data +n;//E(n)
			*x=n-dims[0]+1;
		}
		return PyArray_Return(e);
	}
	PyObject* PyETraj2Hist(PyObject* _traj,PyObject* _Ehist)
	{
		PyArrayObject* traj = (PyArrayObject*) PyArray_ContiguousFromObject(_traj, NPY_INT, 0, 0);
		NULL_CHECK(traj,"PyTraj2Hist(traj:numpy.int)");
		NDIM_CHECK(traj,1);
		
		PyArrayObject* Ehist = (PyArrayObject*) PyArray_ContiguousFromObject(_Ehist, NPY_INT, 0, 0);
		NULL_CHECK(Ehist,"PyTraj2Hist(Ehist:numpy.int)");
		DIM_CHECK(Ehist,0,self->nEn);
			            
        int *trajdata =(int*)traj->data;
        int *Ehistdata =(int*)Ehist->data;
        
        for (int n=0;n<self->nEn;n++)
        	Ehistdata[n]=0;
        	
        for (int c=0; c<traj->dimensions[0]; c++)
            Ehistdata[self->nEn - 1 + trajdata[c]] += 1;
                 
        int e0 = 0, e1 = 0;
        for (int c=0;c<self->nEn;c++)
        {
                if (Ehistdata[c] > 0)
                {
                    e0 = c - self->nEn + 1;
                    break;
                 }
        }
        for (int c=0;c<self->nEn;c++)
        {
                if (Ehistdata[self->nEn - c - 1] > 0)
                {
 					e1 = -c;
                    break;
                 }
        }
        return Py_BuildValue("(iiO)",e0,e1,Ehist); 
	}
	PyObject* PySetWolff()
	{
		PyArrayObject *stack; npy_intp dims[1];
		dims[0]  = self->nN;
		stack = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
		self->ClusterStack_p=(int *)stack->data;
		self->dPadd=1-exp(-1.0/self->Arg_io[self->AI_TEMPER]);
		return PyArray_Return(stack);
	}
	PyObject* PySetGCEWolff()
	{
		PyArrayObject *Seq_p,*Stack_p; npy_intp dims[1];
		dims[0]  = self->nN;
		Stack_p = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
		self->ClusterStack_p=(int *)Stack_p->data;
		Seq_p = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
		self->ClusterSeq_p=(int *)Seq_p->data;
		return Py_BuildValue("(OO)",Stack_p,Seq_p);
	}
	PyObject* PySetWl(double lnfPrecision, double flatRate)
	{
		PyArrayObject *lnge,*gen;
		npy_intp dims[1];
		dims[0] = self->nN+self->nN+1;
		lnge = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		gen = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_ULONG);
		self->SetWL((double *)lnge->data,(unsigned long *)gen->data, gen->dimensions[0], lnfPrecision,  flatRate);
		return Py_BuildValue("(OO)",lnge,gen);
	}
	PyObject* PySetMul(double h)
	{
		PyArrayObject *beta,*se,*gn,*hist;
		npy_intp dims[1];
		dims[0] = self->nN+self->nN+1;
		beta = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		se = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		gn = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		hist = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_ULONG);
		self->SetMul((double *)beta->data,(double *)se->data,(double *)gn->data,(unsigned long *)hist->data, beta->dimensions[0], h);
		return Py_BuildValue("(OOOO)",beta,se,gn,hist);
	}
	PyObject* PyMulReweight(PyObject* _traj)
	{
		PyArrayObject* traj = (PyArrayObject*) PyArray_ContiguousFromObject(_traj, NPY_INT, 0, 0);
		NULL_CHECK(traj,"PyMulReweight(traj:numpy.int)");
		NDIM_CHECK(traj,1);
		int *trajdata =(int*)traj->data;
		for (int c=0; c<traj->dimensions[0]; c++)
            self->Hist_p[self->nEn - 1 + trajdata[c]] += 1;
		if (self->MulReweight())
			return Py_True;
		else
			return Py_False;
	}
	/*
	PyObject* PyReSetFun(PyObject* a_)
  	{
  		PyArrayObject* a = (PyArrayObject*) PyArray_ContiguousFromObject(a_, NPY_DOUBLE, 0, 0);
		NULL_CHECK(a,"PyReSetFun:PyArray_ContiguousFromObject(numpy.double)");
    		NDIM_CHECK(a,2);
    		DIM_GT_CHECK(a,1,5);
    		self->ReSetFun(a->dimensions[0],a->dimensions[1],(double *)a->data);
    		return PyArray_Return(a);
  	}
  	PyObject* PyCalBS(double lflatex,double hflatex, double high_t,double low_t,bool smooth)
	{
		PyArrayObject *B,*S;
		npy_intp dims[1];
		dims[0] = self->nN+self->nN+1;
		B = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		NULL_CHECK(B,"PyCalBS:PyArray_SimpleNew(numpy.double)");
		S = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		NULL_CHECK(S,"PyCalBS:PyArray_SimpleNew(numpy.double)");
		double *f,beta1;
		self->InitFun(lflatex,hflatex,smooth);
		for(int n=0;n<dims[0];n++)
		{
			f=(double *)S->data +n;
			*f=self->SEFunction(n-dims[0]+1,&beta1,high_t,low_t);
			f=(double *)B->data +n;//B(e)
			*f=beta1;
		}
		self->SetSefun((double *)S->data,S->dimensions[0]);
		return Py_BuildValue("(OO)",B,S);
	}
	*/
	PyObject* PySetSe(PyObject *S_)
	{
	  		PyArrayObject* S = (PyArrayObject*) PyArray_ContiguousFromObject(S_, NPY_DOUBLE, 0, 0);
    		self->SetSefun((double *)S->data,S->dimensions[0]);
    		return PyArray_Return(S);
	}
}

%extend SParticle
{
	PyObject* PyGetXY_Pak(double step)
	{
		PyArrayObject *X,*Y; npy_intp dims[1];
		dims[0]  = int((self->PAK_XH-self->PAK_XL)/step);

		X = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		Y = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		double  *x,*y;
		
		double old=self->dX;
		double oldE=self->Arg_io[self->AI_HAMI];
		for(int n=0;n<dims[0];n++)
		{
			self->dX=self->PAK_XL+n*step;
			self->Gauge_Pak();
			x=(double *)X->data +n;
			y=(double *)Y->data +n;
			*y=self->Arg_io[self->AI_HAMI];
			*x= self->dX;
		}
		self->dX=old;
		self->Arg_io[self->AI_HAMI]=oldE;
		return Py_BuildValue("(OO)",X,Y);
	}
	PyObject* PyGetXY_MP(double step)
	{
		PyArrayObject *X,*Y; npy_intp dims[1];
		dims[0]  = int((self->MP_XH-self->MP_XL)/step);

		X = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		Y = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		double  *x,*y;
		
		double old=self->dX;
		double oldE=self->Arg_io[self->AI_HAMI];
		for(int n=0;n<dims[0];n++)
		{
			self->dX=self->MP_XL+n*step;
			self->Gauge_MP();
			x=(double *)X->data +n;
			y=(double *)Y->data +n;
			*y=self->Arg_io[self->AI_HAMI];
			*x= self->dX;
		}
		self->dX=old;
		self->Arg_io[self->AI_HAMI]=oldE;
		return Py_BuildValue("(OO)",X,Y);
	}

}     	
%extend Ising
{
	PyObject* PySetConf(PyObject* a_)
	{
		PyArrayObject* a = (PyArrayObject*) PyArray_ContiguousFromObject(a_, NPY_INT8, 0, 0);
		NULL_CHECK(a,"PySetConf:PyArray_ContiguousFromObject(numpy.int8)");
		NDIM_CHECK(a,2);
		self->SetConf(a->data,a->dimensions[0],a->dimensions[1]);
		return PyArray_Return(a);
	}
	PyObject* PySetWolff()
	{
		PyArrayObject *stack; npy_intp dims[1];
		dims[0]  = self->nN;
		stack = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
		self->ClusterStack_p=(int *)stack->data;
		self->dPadd=1-exp(-2.0/self->Arg_io[self->AI_TEMPER]);
		return PyArray_Return(stack);
	}
	PyObject* PySetGCEWolff()
	{
		PyArrayObject *Seq_p,*Stack_p; npy_intp dims[1];
		dims[0]  = self->nN;
		Stack_p = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
		self->ClusterStack_p=(int *)Stack_p->data;
		Seq_p = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
		self->ClusterSeq_p=(int *)Seq_p->data;
		return Py_BuildValue("(OO)",Stack_p,Seq_p);
	}
	PyObject* PyGetEnergyRange()
	{
		PyArrayObject *e; npy_intp dims[1];
		dims[0]  = self->nN+1;
		e = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_INT);
		int *x;
		for(int n=0;n<dims[0];n++)
		{
			x=(int *)e->data +n;//E(n)
			*x= 4*n-2*self->nN;
		}
		return PyArray_Return(e);
	}
	PyObject* PySetWl(double lnfPrecision, double flatRate)
	{
		PyArrayObject *lnge,*gen;
		npy_intp dims[1];
		dims[0] = self->nN+1;
		lnge = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		gen = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_ULONG);
		self->SetWL((double *)lnge->data,(unsigned long *)gen->data, gen->dimensions[0], lnfPrecision,  flatRate);
		return Py_BuildValue("(OO)",lnge,gen);
	}
	 // hisE[(2N + hami) / 4)]++
	PyObject* PyETraj2Hist(PyObject* _traj,PyObject* _hist)
	{
		PyArrayObject* traj = (PyArrayObject*) PyArray_ContiguousFromObject(_traj, NPY_INT, 0, 0);
		NULL_CHECK(traj,"PyETraj2Hist(traj:numpy.int)");
		NDIM_CHECK(traj,1);
		
		PyArrayObject* hist = (PyArrayObject*) PyArray_ContiguousFromObject(_hist, NPY_INT, 0, 0);
		NULL_CHECK(hist,"PyETraj2Hist(hist:numpy.int)");
		DIM_CHECK(hist,0,self->nEn);
			            
        int *trajdata =(int*)traj->data;
        int *histdata =(int*)hist->data;
        
        for (int n=0;n<self->nEn;n++)
        	histdata[n]=0;
        	
        for (int c=0; c<traj->dimensions[0]; c++)
            histdata[(2*self->nN  + trajdata[c])/4] ++;
                 
        int e0 = 0, e1 = 0;
        for (int c=0;c<self->nEn;c++)
        {
                if (histdata[c] > 0)
                {
                    e0 = 4*c -2* self->nN;
                    break;
                 }
        }
        for (int c=0;c<self->nEn;c++)
        {
                if (histdata[self->nEn - c - 1] > 0)
                {
 					e1 =2*self->nN -4*c;
                    break;
                 }
        }
        return Py_BuildValue("(iiO)",e0,e1,hist); 
	}
	PyObject* PyMTraj2Hist(PyObject* _traj,PyObject* _hist)
	{
		PyArrayObject* traj = (PyArrayObject*) PyArray_ContiguousFromObject(_traj, NPY_INT, 0, 0);
		NULL_CHECK(traj,"PyTraj2Hist(traj:numpy.int)");
		NDIM_CHECK(traj,1);
		
		int len=self->nN+1;
		PyArrayObject* hist = (PyArrayObject*) PyArray_ContiguousFromObject(_hist, NPY_INT, 0, 0);
		NULL_CHECK(hist,"PyMTraj2Hist(hist:numpy.int)");
		DIM_CHECK(hist,0,len);
			            
        int *trajdata =(int*)traj->data;
        int *histdata =(int*)hist->data;

        for (int n=0;n<len;n++)
        	histdata[n]=0;
        	
        for (int c=0; c<traj->dimensions[0]; c++)
            histdata[(self->nN  + trajdata[c])/2] ++;
                 
        int e0 = 0, e1 = 0;
        for (int c=0;c<len;c++)
        {
                if (histdata[c] > 0)
                {
                    e0 = 2*c- self->nN;
                    break;
                 }
        }
        for (int c=0;c<len;c++)
        {
                if (histdata[len - c-1] > 0)
                {
 					e1 =self->nN -2*c;
                    break;
                 }
        }
        return Py_BuildValue("(iiO)",e0,e1,hist); 
	}

}
%extend LqFluid
{
	PyObject* PyGetCor(int ndim0,int ndim1)
	{
		PyArrayObject *Corxs,*Corys,*Fxs, *Fys;
		npy_intp dims[1];
		dims[0]  = ndim0*ndim1;
		Corxs = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		Corys = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		Fxs = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		Fys = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		NULL_CHECK(Corxs,"PyGetCor,Corxs");
		NULL_CHECK(Corys,"PyGetCor,Corys");
		NULL_CHECK(Fxs,"PyGetCor,Fxs");
		NULL_CHECK(Fys,"PyGetCor,Fys");
		self->SetConf((double *)Corxs->data, (double *)Corys->data, (double *)Fxs->data, (double *)Fys->data, ndim0,ndim1);
    	return Py_BuildValue("(OOOO)",Corxs,Corys,Fxs, Fys);
	}
	PyObject* PyGetXY_LQ(double begin,double end,double num)
	{
		PyArrayObject *X,*Y;
		npy_intp dims[1];
		dims[0]  =num;
		double dx=(end-begin)/num;

		X = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		Y = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		double  *x,*y;

		for(int n=0;n<dims[0];n++)
		{
			x=(double *)X->data +n;
			y=(double *)Y->data +n;
			*x=begin+n*dx;
			*y=self->GetLqValue(*x);
		}
		return Py_BuildValue("(OO)",X,Y);
	}
	
	PyObject* PyGetRDF(PyObject* a_)
	{
		PyArrayObject* a = (PyArrayObject*) PyArray_ContiguousFromObject(a_, NPY_DOUBLE, 0, 0);
		NULL_CHECK(a,"PyGetRDF(NPY_DOUBLE)");
		int ni=self->GetRDF((double *)a->data,a->dimensions[0]);
		return Py_BuildValue("(Oi)",a,ni);
	}
}

%inline %{
	PyObject* PyAutoCor(PyObject* a_,int kmax)
  	{
  		PyArrayObject* a = (PyArrayObject*) PyArray_ContiguousFromObject(a_, NPY_DOUBLE, 0, 0);
		if (a==NULL || a->nd !=1) {
			PyErr_Format(PyExc_TypeError, 	"%s is NULL, %s", "PyArrayObject* a","PyAutoCor"); 
          	return NULL;
     	}

		PyArrayObject *Xt; npy_intp dims[1];
		dims[0]  = kmax+1;
		Xt = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  		int ta=AutoCor((double *)a->data, a->dimensions[0], (double *)Xt->data, kmax);		
    	return Py_BuildValue("(iO)",ta,Xt);
  	}
  	
  	PyObject*  PyMeanStdErr(PyObject* TR,  unsigned int tau) 
  	{
		PyArrayObject* a = (PyArrayObject*) PyArray_ContiguousFromObject(TR, NPY_DOUBLE, 0, 0);
		if (a==NULL || a->nd !=1) {
			PyErr_Format(PyExc_TypeError, 	"%s is NULL, %s", "PyArrayObject* a","PyMeanStdErr"); 
          	return NULL;
     	}
     	double mean,std,meanErr; 	  		
		MeanStdErr((double *)a->data, a->dimensions[0], &mean, &std,&meanErr,tau);
		
		return Py_BuildValue("(ddd)",mean,std,meanErr);
	}
  	PyObject*  PyJackknife_MeanStdErrTau_I(PyObject* TR,  unsigned int nBlock,  unsigned int MaxTau) 
  	{
		PyArrayObject* a = (PyArrayObject*) PyArray_ContiguousFromObject(TR, NPY_INT, 0, 0);
		if (a==NULL || a->nd !=1) {
			PyErr_Format(PyExc_TypeError, 	"%s is NULL, %s", "PyArrayObject* a","PyJackknife_MeanStdErrTau_I"); 
          	return NULL;
     	}
     	PyArrayObject *tauCorr; npy_intp dims[1];
		dims[0]  = MaxTau;
		tauCorr = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		
     	double mean=0,std=0,meanErr=0,ptau=0;
		int taui=Jackknife_MeanStdErrTau_I((int *)a->data,a->dimensions[0],nBlock,(double *)tauCorr->data, MaxTau, &mean, &std,&meanErr,&ptau);
		
		return Py_BuildValue("(ddddiO)",mean,std,meanErr,ptau,taui,tauCorr);
	}
	PyObject*  PyJackknife_MeanStdErrTau_D(PyObject* TR,  unsigned int nBlock,  unsigned int MaxTau) 
  	{
		PyArrayObject* a = (PyArrayObject*) PyArray_ContiguousFromObject(TR, NPY_DOUBLE, 0, 0);
		if (a==NULL || a->nd !=1) {
			PyErr_Format(PyExc_TypeError, 	"%s is NULL,%s", "PyArrayObject* a", "PyJackknife_MeanStdErrTau_D"); 
          	return NULL;
     	}
     	PyArrayObject *tauCorr; npy_intp dims[1];
		dims[0]  = MaxTau;
		tauCorr = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
		
     	double mean=0,std=0,meanErr=0,ptau=0;
		int taui=Jackknife_MeanStdErrTau_D((double *)a->data,a->dimensions[0],nBlock,(double *)tauCorr->data, MaxTau, &mean, &std,&meanErr,&ptau);
		
		return Py_BuildValue("(ddddiO)",mean,std,meanErr,ptau,taui,tauCorr);
	}
%}

%init%{
    import_array();
%}
