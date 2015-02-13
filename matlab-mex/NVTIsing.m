function [out]=NVTIsing(n,m,temp,samples)
conf=zeros(n,m,'int8');
ising('new',conf);
ising('InitConf',0);
ising('SetTemperature',temp);
%E=ising('GetEnergyRange');
TR=zeros(1,samples);
stride=n*m*10;
tic;
for i=1:samples
    TR(i)=ising('NVT',stride);
end
toc
disp(ising('GetNVTArg_io'));

PlotStat(TR,n,m);

out(1)=temp;
out(2)=mean(TR);%E
out(3)=std(TR);%sigma

%clear ising
end

function PlotStat(TR,n,m)
%figure;
subplot(211);
hist(TR,n*m*2);
xlabel('E');
ylabel('Visited Numbers');
grid on;

subplot(212);
plot(TR);
xlabel('Samples');
ylabel('Trajectories of E');
grid on;
end

