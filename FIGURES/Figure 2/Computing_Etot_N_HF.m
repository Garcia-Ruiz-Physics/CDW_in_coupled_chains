%% Total energy for two chains: NUMERICAL APPROACH
% This programme computes the total energy for 
% two parallel-coupled chains with a CDW order.

%% clear the workspace
clear all
close all
clc
%% Defining constants
t=1;
a=1;
%% Defining problem parameters
stepD=10^-4;
MaxD=0.1;
D=0:stepD:MaxD;
lD=length(D);
stepT=0.2*10^-3;
tperp=0:stepT:0.025;
lT=length(tperp);
for N=6 % Number of chains
    E_electronic=zeros(lT,lD);
    for T=1:lT
        %% COUPLING-INDUCED SHIFTS
        SHIFT=zeros(1,N);
        for n=1:N
            SHIFT(n)=2*tperp(T)*cos(n*pi/(N+1));
        end
        %% Get the job done
        % momentum grid
        stepK=2*10^-6;
        k=(-1:stepK:0)*2*pi/(2*a)/2;
        lK=length(k);
        % At half filling, the Fermi level lies always at zero!
        % Band structure
        Ee=zeros(1,lD);
        for DDD=1:lD
            Ek=2*t*sqrt(1-(1-D(DDD)^2)*( sin(k*a) ).^2);
            Ek_bands=zeros(N,lK);
            for n=1:N
                Ek_bands(n,:)=Ek+SHIFT(n);
            end
            Em_bands=-Ek_bands;
            %% Computing the total energy
            Ee(DDD)=2*...% We have only computed half of the BZ
                sum( ... % Sum over all SHIFTS
                sum(Ek_bands.*(Ek_bands<0))+ ... % Positive bands
                sum(Em_bands.*(Em_bands<0))  ... % Negative bands
                )/(2*lK); % Total number of unit cells
        end
        E_electronic(T,:)=Ee;
        disp([num2str(T/lT*100),' % done'])
    end
    save(['Electronic_contribution_Half_Filled_N_',num2str(N),'_v2'], ...
         'N','E_electronic','tperp','D','stepK')
end
