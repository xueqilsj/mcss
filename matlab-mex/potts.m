%POTTS Monte Carlo Statistical Simulation based on MATLAB+MEX.
%   POTTS for NVT,Wang-Landau,Fast Adaptive Flat Histogram.
%
%	POTTS(conf,Q) ,conf=zeros(n,m,'int8'),Q>0
%	POTTS("SetTemperature",t),t>0
%	POTTS("InitConf",d),d=0:random,d<0:ground state,d>0:anti-ground state
%	POTTS("Gauge")
%	[Arg]=POTTS("GetNVTArg_io")
%
%	[HAMI]=POTTS("NVT",trials);trials>0
%
%	[HAMI]=POTTS( "WL",trials)
%	[LnGe, GeN]=POTTS('SetWL',lnfPrecision, flatRate);lnfPrecision=1e-8, flatRate=0.90
%	POTTS("WLUpdatePrecision",p);p<1
%	E=POTTS("GetEnergyRange")
%	POTTS("WLCheckBelowAVG")
%	POTTS("WLResetGeN")
%
%	[HAMI]=POTTS("FAFH",trials)
%
%   See also ...
%   MEX-File function.


