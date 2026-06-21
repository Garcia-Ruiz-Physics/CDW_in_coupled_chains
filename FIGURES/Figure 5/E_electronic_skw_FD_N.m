%% Total energy for multiple skewed-coupled chains: NUMERICAL APPROACH
% This programme computes the total energy for 
% two skewed-coupled chains with a CDW order.

%% clear the workspace
clear all
close all
clc
%% Defining constants
t=1;
a=1;
for N=3:6
    %% Defining problem parameters
    stepD=5*10^-4;
    % Peierls CDW phase
    D=[];
    D0=0.1:0.1:9;
    for OoM=-4:2:-2
        D=[D D0*10^OoM];
    end
    D=[0 D 0.091:0.001:0.1];
    lD=length(D);
    % Inter-chain coupling
    tperp=0.02;
    lT=length(tperp);
    % momentum grid
    stepK=5*10^-6;
    % stepK=50*10^-6; % YOU CAN CHANGE IT FOR BETTER ACCURACY
    k=(-1:stepK:0)*2*pi/(2*a)/2;
    lK=length(k);
    % Doping levels
    doping=0:0.00002:0.01;
    % doping=0:0.00005:0.01; % YOU CAN CHANGE IT FOR BETTER ACCURACY
    lDO=length(doping);
    Ntotal=lK*N*2;
    doping2points=doping... Proportional to doping level
        /0.5* ... would be equivalent to complete filling
        (lK*N)+... The remaining number of points
        lK*N; % If three lines above are zero, we are at half filling

    doping2points=round(doping2points);
    Ee_nem=zeros(lT,lD,lDO);
    Ee_zzg=zeros(lT,lD,lDO);
    for T=1:lT % For every inter-chain coupling
        for d=1:lD % For every Peierls CDW order
            Enem=zeros(lK,2*N);
            Ezzg=zeros(lK,2*N);
            tic
            for K=1:lK
                % For every momentum point
                %% 1D atomic chain
                f=(t+D(d))+(t-D(d))*exp(2i*k(K)*a);
                H_1D=[0  f;
                    f' 0];
                Vs=[1 exp(2i*k(K)*a);
                    1 1                   ]*tperp(T);
                %% Decoupled Hamiltonian
                H_dec=[];
                for n=1:N
                    H_dec=blkdiag(H_dec,H_1D);
                end
                %% Nematic phase
                T_nem=[];
                for NNN=1:N-1
                    T_nem=blkdiag(T_nem,Vs);
                end
                T_nem=[zeros(2*(N-1),2) T_nem;
                    zeros(2,2*N)];
                T_nem=T_nem+T_nem';
                H_nem=H_dec+T_nem;
                %% Zig-Zag phase
                T_zzg=[];
                for NNN=1:N-1
                    if NNN/2==floor(NNN/2)
                        T_zzg=blkdiag(T_zzg,Vs');
                    else
                        T_zzg=blkdiag(T_zzg,Vs);
                    end
                end
                T_zzg=[zeros(2*(N-1),2) T_zzg;
                    zeros(2,2*N)];
                T_zzg=T_zzg+T_zzg';
                H_zzg=H_dec+T_zzg;
                %% Diagonalizing the Hamiltonian
                Enem(K,:)=eig(H_nem);
                Ezzg(K,:)=eig(H_zzg);
            end
            Values_nem=reshape(Enem,1,lK*2*N);
            Values_zzg=reshape(Ezzg,1,lK*2*N); %% THERE WAS A TYPO HERE...
            toc
            %% Computing the EF for nem
            tic
            for DO=1:lDO
                [Energy_values_sorted,~]=sort(Values_nem);
                EF_nem=(Energy_values_sorted( doping2points(DO))+...
                    Energy_values_sorted( doping2points(DO)+1) )/2;
                [Energy_values_sorted,~]=sort(Values_zzg);
                EF_zzg=(Energy_values_sorted( doping2points(DO))+...
                    Energy_values_sorted( doping2points(DO)+1) )/2;
                %% Total electronic energy up to some kc value (see k)
                Total_Enem_list=sum( Values_nem.*(Values_nem<EF_nem) );
                Total_Ezzg_list=sum( Values_zzg.*(Values_zzg<EF_zzg) );%% THERE WAS A TYPO HERE...
                Ee_nem(T,d,DO)=2*...% Two "valleys"
                    (Total_Enem_list)... % Sum over all bands
                    /(2*lK); % Total number of unit cells
                Ee_zzg(T,d,DO)=2*...% Two "valleys"
                    (Total_Ezzg_list)... % Sum over all bands %% THERE WAS A TYPO HERE...
                    /(2*lK); % Total number of unit cells
            end
            toc
            disp([num2str(d/lD*100),' % done'])
        end
        % disp([num2str(T/lT*100),' % done'])
    end
    save(['Electronic_contribution_Finite_Doping_skewed_N', ...
        num2str(N),'_tperp_',num2str(1000*tperp),'e-3t'], ...
        'N','Ee_nem','Ee_zzg','tperp','D','stepK','doping')
end