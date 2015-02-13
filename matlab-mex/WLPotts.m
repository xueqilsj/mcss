function WLPotts(n,m,q,lnfPrecision,flatRate)%1e-8, flatRate=0.90
conf=zeros(n,m,'int8');
potts('new',conf,q);
potts('InitConf',0);
potts('Gauge');
E = potts('GetEnergyRange');
[LnGe, GeN] = potts('SetWl',lnfPrecision, flatRate);
stride=n*m*12;
tic;
tstart=tic;
samples = 0;
potts('WLUpdatePrecision',1);
Lnf=1;
%format long g
while Lnf > lnfPrecision
    samples = samples+1;
    Tr(samples)=potts('WL',stride);
    if potts('WLCheckBelowAVG')
        Lnf = Lnf/2;
        fprintf('%d\t%8.4f\t%g\n',samples,toc(tstart),potts('WLUpdatePrecision',Lnf));
        if Lnf > lnfPrecision
            samples = 0;
            clear Tr;
            potts('WLResetGeN');
        end
        
    end
end
toc;

fprintf('mean(GeN)=%f\n',mean(GeN));
PlotStat(Tr,n,m,q,ExactRange(E),ExactRange(LnGe),ExactRange(GeN));

%clear potts
end

function PlotStat(Tr,n,m,q,E,LnGe,GeN)
subplot(211);
plot(E,GeN);
xlabel('E');
ylabel('Number of G(E)');
title(sprintf('WLPotts%dx%d,Q=%d',n,m,q));
grid on;

subplot(212);
hist(Tr,n*m*2);
xlabel('Samples');
ylabel('Trajectories of E');
grid on;

figure;
plot(E,LnGe);
xlabel('E');
ylabel('LnG(e)');
title(sprintf('WLPotts%dx%d,Q=%d',n,m,q));
grid on;

end

function [out]=ExactRange(in)
out = in(5:end); %trim 4 unaccessible energy levels
out(1) = in(1);
out(2) = in(5);
end

