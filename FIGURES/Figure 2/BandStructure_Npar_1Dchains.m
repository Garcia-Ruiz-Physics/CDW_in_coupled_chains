%% Total energy for two chains: NUMERICAL APPROACH
% This programme computes the total energy for 
% two parallel-coupled chains with a CDW order.

function BandStructure_Npar_1Dchains(N,tperp)
%% Defining constants
t=1;
a=1;
%% Defining problem parameters
MaxD=0.004;
D=[ 0 MaxD];
%% COUPLING-INDUCED SHIFTS
SHIFT=zeros(1,N);
for n=1:N
    SHIFT(n)=2*tperp*cos(n*pi/(N+1));
end
%% Get the job done
% momentum grid
stepK=2*10^-6;
k=(-1:stepK:0.9)*2*pi/(2*a)/2;
% At half filling, the Fermi level lies always at zero!
% Band structure
hold off
Ek_CDW=2*t*sqrt(1-(1-D(1)^2)*( sin(k*a) ).^2);
for n=1:N
    Ek_CDW_bands=Ek_CDW+SHIFT(n);
    plot(k/(pi),100*Ek_CDW_bands,'b','LineWidth',1)
    hold on
    plot(k/(pi),-100*Ek_CDW_bands,'b','LineWidth',1)
end


Ek_CDW=2*t*sqrt(1-(1-D(2)^2)*( sin(k*a) ).^2);
for n=1:N
    Ek_CDW_bands=Ek_CDW+SHIFT(n);
    plot(k/(pi),100*Ek_CDW_bands,'-','color',0.8*[1 1 0],'LineWidth',1.5)
    hold on
    plot(k/(pi),-100*Ek_CDW_bands,'-','color',0.8*[1 1 0],'LineWidth',1.5)
end

Ek_CDW=2*t*sqrt(1-(1-D(1)^2)*( sin(k*a) ).^2);
for n=1:N
    Ek_CDW_bands=Ek_CDW+SHIFT(n);
    plot(k/(pi),100*Ek_CDW_bands,'b--','LineWidth',1)
    hold on
    plot(k/(pi),-100*Ek_CDW_bands,'b--','LineWidth',1)
end
xlim([-0.5 -0.49])
ylim([-1 1]*0.1)
