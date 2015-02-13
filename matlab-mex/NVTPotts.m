function [out]=NVTPotts(n,m,q,temp,samples)
conf=zeros(n,m,'int8');
potts('new',conf,q);
potts('InitConf',0);
potts('SetTemperature',temp);
%E=potts('GetEnergyRange');
TR=zeros(1,samples);
stride=n*m*10;
tic;
for i=1:samples
    TR(i)=potts('NVT',stride);
end
toc
disp(potts('GetArg_io'));

PlotStat(TR,n,m);

out(1)=temp;
out(2)=mean(TR);%E
out(3)=std(TR);%sigma

%clear potts
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
ylabel('TRajectories of E');
grid on;
end

