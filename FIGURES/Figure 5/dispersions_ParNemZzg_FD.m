%% Band structure of N skewed-coupled atomic chains: numerical
function dispersions_ParNemZzg_FD(N,tperp)
%% INPUTS:
% N: number of chains
% tperp: Inter-chain coupling
%% Defining constants
green =[19 190 184]/255;
% yellow=[240 240 20]/255;
yellow=[1 1 0]*0.8;
e=1; % Electron's charge
a=1;
t=e; % Hopping parameter
%% CDW order parameters
aX=[0 5]*10^-3*t; % intra-chain hopping parameter difference (2 points version)
laX=length(aX);
%% INTER-CHAIN COUPLING parameters
tperp=tperp*t; % We should fix this parameter at 0.025t, but it could be an INPUT
Vp=eye(2)*tperp;
%% Defining the momentum mesh: First choice (constant)
stepK=10^-5; % Change for better accuracy!
minK=-0.5;
maxK=-0.49; % At most you should consider 0
k_list=(minK:stepK:maxK)*2*pi/(2*a);
lK=length(k_list);
%% Vectorising shitte
Epar=zeros(laX,lK,2*N);
Enem=zeros(laX,lK,2*N);
Ezzg=zeros(laX,lK,2*N);
for A=1:laX
    for KKK=1:lK
        %% 1D atomic chain
        f=(t+aX(A)/2)+(t-aX(A)/2)*exp(2i*k_list(KKK)*a);
        H_1D=[0  f;
            f' 0];
        Vs=[1 exp(2i*k_list(KKK)*a);
            1 1               ]*tperp;
        %% Decoupled Hamiltonian
        H_dec=[];
        for n=1:N
            H_dec=blkdiag(H_dec,H_1D);
        end
        %% Parallel phase
        T_par=[];
        for NNN=1:N-1
            T_par=blkdiag(T_par,Vp);
        end
        T_par=[zeros(2*(N-1),2) T_par;
            zeros(2,2*N)];
        T_par=T_par+T_par';
        H_par=H_dec+T_par;
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
        %% Diagonalizing shitte
        Epar(A,KKK,:)=eig(H_par);
        Enem(A,KKK,:)=eig(H_nem);
        Ezzg(A,KKK,:)=eig(H_zzg);
    end
    %% Graphical output: band structure
    if A==1
        hold on
        for B=1:2*N
            plot(k_list/(pi),Enem(1,:,B),'b-','LineWidth',0.5)
        end
        xlim([min(k_list) max(k_list)]/pi)
        ylim([-1 1]*0.06)
    elseif A==laX
        hold on
        for B=1:2*N
            plot(k_list/pi,Enem(A,:,B),'color',green,'LineWidth',2.5)
            plot(k_list/pi,Ezzg(A,:,B),'color',yellow,'LineWidth',2)
        end
        box on
        xlim([min(k_list) max(k_list)]/pi)
        ylim([-1 1]*0.06)
        xticks([-0.5 -0.49])
        yticks([-0.05 0 0.05])

        for B=1:2*N
            plot(k_list/(pi),Enem(1,:,B),'b--','LineWidth',0.5)
        end
    end
end
end