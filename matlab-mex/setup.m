function setup(varargin)
%COREDIR = ['..' filesep 'core-cpp' filesep];
if nargin>0 && strcmp(varargin{1},'clean')==1
    delete('*.mexglx');
    delete('*.mexw32');
else
    mex -I../core-cpp/ -L../core-cpp -output potts Potts_mex.cpp ../core-cpp/Potts.cpp ../core-cpp/SeFun.cpp ../core-cpp/Random.cpp
    mex -I../core-cpp/ -L../core-cpp -output mclinear Linear_mex.cpp ../core-cpp/Linear.cpp ../core-cpp/SeFun.cpp ../core-cpp/Random.cpp
    mex -I../core-cpp/ -L../core-cpp -output ising Ising_mex.cpp ../core-cpp/Ising.cpp ../core-cpp/SeFun.cpp ../core-cpp/Random.cpp
end
end