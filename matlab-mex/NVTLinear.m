function [out]=NVTLinear(initp, offset, temp,samples)
mclinear('InitParameters',initp, offset, temp);
[X,Y]=mclinear('GetXY',0.01);
TR=zeros(1,samples);
TRX=zeros(1,samples);
stride=5000;
tic;
for i=1:samples
    TR(i)=mclinear('NVT',stride);
    TRX(i)=mclinear('GetPoint',0);
end
toc
disp(mclinear('GetArg_io'));

PlotStat(X,Y,TR,TRX);

out(1)=temp;
out(2)=mean(TR);%E
out(3)=std(TR);%sigma

%clear mclinear
end

function PlotStat(X,Y,TR,TRX)
%figure;
subplot(311);
plot(X, Y);
xlabel('x');
ylabel('F(x)');
grid on;

subplot(312);
%x=-2.0:0.004:2.0;n = hist(TRX, x);
plot(TRX);
xlabel('Samples');
ylabel('Histogram');
grid on;

subplot(313);
plot(TR);
xlabel('Samples');
ylabel('Traj. of E');
grid on;
end

