%% Figure 3
% This programme produces Fig. 3 of the paper:
% "Dimensional and Doping Stability of Peierls CDW
%  in Arrays of Coupled One-Dimensional Chains"

%% Clean the workspace
clear all
close all
clc
% Defining constants
X0=[0.565 0.645 0.69 0.72];
Y0=[0.858 0.804 0.761 0.717];
%% Figure positioning
f=figure;
f.Position=[150,150,[1500,900]];
for N=3:6
    kappa=0:0.001:1.5;
    lKA=length(kappa);
    %% Load data
    load(['Electronic_contribution_Half_Filled_N_',...
        num2str(N),'_v2.mat'])
    lT=length(tperp);
    Epot=4*N*D'.^2*kappa; % Elastic energy
    % For different values of tperp, you have different maps
    Minima=zeros(lT,lKA);
    for Nt=1:lT
        flag=0;
        for K=1:lKA
            Etot=(E_electronic(Nt,:)+(Epot(:,K))'+4/pi*N)*10^3;
            [nMin, idxMin, valMin] = localMinima(Etot);
            if nMin==1 && idxMin>1
                Minima(Nt,K)=2;
            elseif nMin>1
                flag=1;
                Minima(Nt,K)=1;
            end
            if nMin==1 && flag==1
                if N/2==floor(N/2)
                    Minima(Nt,K)=0;
                else
                    Minima(Nt,K)=0.5;
                end
            end
        end
    end

    %% %%%%%%%%%%%% %%
    %% UPPER PANELS %%
    %% %%%%%%%%%%%% %%
    axes('Position',[ 0.06+0.22*(N-3) 0.58 0.20 0.4])
    surf(tperp,kappa,Minima','EdgeColor','none')
    hold on
    view([0 90])
    if N==3
        ylabel('\kappa (t^{-1})')
    else
        yticks([])
    end
    xlabel('$\gamma$ ($t$)','Interpreter','latex')
    %% Analytical solution
    % if N/2==floor(N/2)
    kappa_an=1/(2*pi)*(log(cos(tperp/4)./sin(tperp/4))-cos(tperp/2));
    plot3(tperp,kappa_an,ones(1,lT)*2,'r--','LineWidth',1.5)
    % end
    axis tight
    set(gca,'FontSize',20)

    text(0.019,1.4,2,['N = ',num2str(N)],'Fontsize',20,'Color','w')
    clim([0 2])
    if N==3
        text(-0.007,1.5,'(a)','FontSize',20)
    end
    %% With doping
    kappa=0:0.001:1.5;
    lKA=length(kappa);
    %% Load data
    load(['Electronic_contribution_Finite_Doping_N_', ...
        num2str(N),'_tperp_25e-3t'])
    lDo=length(doping);
    Epot=4*N*D'.^2*kappa; % Elastic energy
    % For different values of tperp, you have different maps
    Minima=zeros(lDo,lKA);
    for Nd=1:lDo
        flag=0;
        for K=1:lKA
            Etot=(E_electronic(Nd,:)+(Epot(:,K))'+4/pi*N)*10^3;
            [nMin, idxMin, valMin] = localMinima(Etot);
            if N/2==floor(N/2)
                if nMin==1 && idxMin>1
                    Minima(Nd,K)=2;
                elseif nMin>1 % Bistability
                    flag=1;
                    Minima(Nd,K)=1;
                end
                if nMin==1 && flag==1 && idxMin==1
                    Minima(Nd,K)=0;
                end
                if nMin==1 && flag==1 && idxMin>1
                    Minima(Nd,K)=0.5;
                end
            else
                if nMin==1 && idxMin>1
                    Minima(Nd,K)=2;
                elseif nMin>1 % Bistability
                    flag=1;
                    Minima(Nd,K)=1;
                end
                if nMin==1 && flag==1 && idxMin==1
                    Minima(Nd,K)=0;
                end
                if nMin==1 && flag==1 && idxMin>1
                    Minima(Nd,K)=0.5;
                end
            end
        end
    end
    %% %%%%%%%%%%%% %%
    %% LOWER PANELS %%
    %% %%%%%%%%%%%% %%
    axes('Position',[ 0.06+0.22*(N-3) 0.08 0.20 0.4])
    surf(100*doping,kappa,Minima','EdgeColor','none')
    hold on
    view([0 90])
    xlabel('n ($\%$)','Interpreter','latex')
    if N==3
        ylabel('\kappa (t^{-1})')
    else
        yticks([])
    end

    if N==3
        text(-0.28,1.5,'(b)','FontSize',20)
    end
    plot3([1 1]*X0(N-2),[0 1]*Y0(N-2),[2 2],'r:','LineWidth',1)
    %% Inset: band structures
    axis tight
    set(gca,'FontSize',20)
    axes('Position',[ 0.18+0.22*(N-3) 0.32 0.07 0.15])
    BandStructure_Npar_1Dchains(N,tperp)

    plot([-0.5 -0.49],100*[1 1]*2*tperp*cos(pi/(N+1)),'r:','LineWidth',1)

    plot([-0.5 -0.49],100*[1 1]*0,'k:','LineWidth',1)
    ylim([-1 1]*6)
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    xlabel('Momentum')
    ylabel('Energy')
    xticks([])
    yticks([])
    set(gca,'FontSize',13)
end
%% A plethora of functions

function [nMin, idxMin, valMin] = localMinima(A)
A = A(:).';          % ensure row vector
N = numel(A);
isMin = false(1,N);
% Interior points
isMin(2:N-1) = A(2:N-1) < A(1:N-2) & A(2:N-1) < A(3:N);
% Endpoints (one-sided)
isMin(1)   = A(1) < A(2);       % left boundary
isMin(N)   = A(N) < A(N-1);     % right boundary
% Outputs
idxMin = find(isMin);
valMin = A(idxMin);
nMin   = numel(idxMin);

end
