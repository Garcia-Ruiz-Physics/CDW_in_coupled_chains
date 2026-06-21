%% E_electronic for N parallel-coupled chains at finite doping
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
tperp=0.025;
for N=3:6 % Number of chains
    %% COUPLING-INDUCED SHIFTS
    SHIFT=zeros(1,N);
    for n=1:N
        SHIFT(n)=2*tperp*cos(n*pi/(N+1));
    end
    %% Get the job done
    % momentum grid
    stepK=2*10^-6;
    k=(-1:stepK:0)*2*pi/(2*a)/2;
    lK=length(k);   
    % Doping levels
    doping=0:0.00005:0.01;
    lDO=length(doping);
    Ntotal=lK*N*2;
    doping2points=doping... Proportional to doping level
                 *(lK*N)... The remaining number of points
                 /0.5 + ... would be equivalent to complete filling
                 lK*N;
    doping2points=round(doping2points);
    % Band structure
    E_electronic=zeros(lDO,lD);
    for DDD=1:lD
        Ek=2*t*sqrt(1-(1-D(DDD)^2)*( sin(k*a) ).^2);
        Ek_bands=zeros(N,lK);
        for n=1:N
            Ek_bands(n,:)=Ek+SHIFT(n);
        end
        Em_bands=-Ek_bands;
        Energy_values=[reshape(Em_bands,1,lK*N) reshape(Ek_bands,1,lK*N)];
        for DO=1:lDO
            [Energy_values_sorted,ORDER]=sort(Energy_values);
            EF=(Energy_values_sorted( doping2points(DO))+...
                Energy_values_sorted( doping2points(DO)+1) )/2;
            %% Computing the total energy
            E_electronic(DO,DDD)=2*...% We have only computed half of the BZ
                sum( ... % Sum over all SHIFTS
                sum(Ek_bands.*(Ek_bands<EF))+ ... % Positive bands
                sum(Em_bands.*(Em_bands<EF))  ... % Negative bands
                )/(2*lK); % Total number of unit cells
        end
        disp([num2str(DDD/lD*100),' % done'])
    end
    save(['Electronic_contribution_Finite_Doping_N_', ...
        num2str(N),'_tperp_',num2str(tperp*1000),'e-3t'], ...
         'N','E_electronic','tperp','doping','D','stepK')
end
