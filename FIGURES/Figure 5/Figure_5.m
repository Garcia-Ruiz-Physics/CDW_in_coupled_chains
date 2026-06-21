%% Figure 5
% This programme produces Fig. 5 of the paper:
% "Dimensional and Doping Stability of Peierls CDW
%  in Arrays of Coupled One-Dimensional Chains"

%% clean the workspace
clear all
close all
clc

%% Set the position of the figure
f=figure;
f.Position=[0,150,[2000,600]];
green =[19 190 184]/255;
yellow=[1 1 0]*0.8;
%% Finite doping
kappa=0:0.001:1.5;
lKA=length(kappa);
for N=3:6
    load(['Electronic_contribution_Finite_Doping_skewed_N', ...
        num2str(N),'_tperp_20e-3t'])
    Ee_nem=squeeze(Ee_nem);
    Ee_zzg=squeeze(Ee_zzg);
    %% Defining shitte
    lDo=length(doping); 
    Epot=4*N*D'.^2*kappa;
    % For different values of tperp, you have different maps
    Minima=zeros(lDo,lKA);
    for Nd=1:lDo
        flag=0;
        for K=1:lKA
            Etot_nem=(Ee_nem(:,Nd)+(Epot(:,K))+4/pi*N)*10^3;
            Etot_zzg=(Ee_zzg(:,Nd)+(Epot(:,K))+4/pi*N)*10^3;
            [nMin_nem, idxMin_nem, valMin_nem] = localMinima(Etot_nem);
            [nMin_zzg, idxMin_zzg, valMin_zzg] = localMinima(Etot_zzg);
            
            if min(valMin_nem)>min(valMin_zzg)
                % IF ZIG-ZAG WINS
                if idxMin_zzg==1
                    Minima(Nd,K)=0;
                else
                    Minima(Nd,K)=2;
                end
            else
                % IF NEMATIC WINS
                if idxMin_nem==1
                    Minima(Nd,K)=0;
                else
                    Minima(Nd,K)=1;
                end
            end
        end
    end
    axes('Position',[ 0.04+0.21*(N-3) 0.12 0.20 0.85])
    surf(100*doping,kappa,Minima','EdgeColor','none')
    hold on
    view([0 90])
    xlabel('n ($\%$)','Interpreter','latex')
    if N==3
        KAPPA=text(-0.08,0.7,'$\kappa (t^{-1})$','Interpreter','latex','FontSize',20,'Color','k');
        set(KAPPA, {'Rotation'}, num2cell(90))
    else
        yticks([])
    end
    set(gca,'FontSize',18)
    text(0.1,1.4,2,['N = ',num2str(N)],'FontSize',20,'Color','w')
    axes('Position',[ 0.17+0.21*(N-3) 0.58 0.06 0.38])
    dispersions_ParNemZzg_FD(N,tperp)
    xlabel('Momentum','Color','w')
    ylabel('Energy','Color','w')
    set(gca,'FontSize',15)

    grid on
    plot([-0.5 -0.49],[0 0],'k:','LineWidth',0.2)
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    xticks([])
    yticks([])
    
end
%% %%%%%%%%%%%%%% %%
%% METALLIC STATE %%
%% %%%%%%%%%%%%%% %%
radius=0.3;
axes('Position',[ 0.88 0.7 0.11 0.3])
for B=0:1
    ball([0 0 2*B],radius,'b')
    hold on
    ball([1 0 2*B],radius,'b')
    ball([2 0 2*B],radius,'b')
    ball([3 0 2*B],radius,'b')
    ball([4 0 2*B],radius,'b')
    ball([5 0 2*B],radius,'b')
end
for B=0:1
    ball([0.5 0 2*B+1],radius,'b')
    hold on
    ball([1.5 0 2*B+1],radius,'b')
    ball([2.5 0 2*B+1],radius,'b')
    ball([3.5 0 2*B+1],radius,'b')
    ball([4.5 0 2*B+1],radius,'b')
    ball([5.5 0 2*B+1],radius,'b')
end
axis tight
axis equal
view([0 0])
camlight(light)
xlim([-0.5 6])
axis off
%% %%%%%%%%%%%%% %%
%% ZIG-ZAG STATE %%
%% %%%%%%%%%%%%% %%
SHIFT=0.13;
radius=0.3;
axes('Position',[ 0.88 0.4 0.11 0.3])
for B=0:3
    ball([0+SHIFT+0.5*mod(B,2) 0 B],radius,yellow)
    hold on
    ball([1-SHIFT+0.5*mod(B,2) 0 B],radius,yellow)
    ball([2+SHIFT+0.5*mod(B,2) 0 B],radius,yellow)
    ball([3-SHIFT+0.5*mod(B,2) 0 B],radius,yellow)
    ball([4+SHIFT+0.5*mod(B,2) 0 B],radius,yellow)
    ball([5-SHIFT+0.5*mod(B,2) 0 B],radius,yellow)
end
axis tight
axis equal
view([0 0])
camlight(light)
xlim([-0.5 6])
axis off
%% %%%%%%%%%%%%% %%
%% NEMATIC STATE %%
%% %%%%%%%%%%%%% %%
SHIFT=0.13;
radius=0.3;
axes('Position',[ 0.88 0.1 0.11 0.3])
    ball([0+SHIFT 0 0],radius,green)
    hold on
    ball([1-SHIFT 0 0],radius,green)
    ball([2+SHIFT 0 0],radius,green)
    ball([3-SHIFT 0 0],radius,green)
    ball([4+SHIFT 0 0],radius,green)
    ball([5-SHIFT 0 0],radius,green)

    ball([0+SHIFT+0.5 0 1],radius,green)
    ball([1-SHIFT+0.5 0 1],radius,green)
    ball([2+SHIFT+0.5 0 1],radius,green)
    ball([3-SHIFT+0.5 0 1],radius,green)
    ball([4+SHIFT+0.5 0 1],radius,green)
    ball([5-SHIFT+0.5 0 1],radius,green)

    ball([0-SHIFT 0 2],radius,green)
    ball([1+SHIFT 0 2],radius,green)
    ball([2-SHIFT 0 2],radius,green)
    ball([3+SHIFT 0 2],radius,green)
    ball([4-SHIFT 0 2],radius,green)
    ball([5+SHIFT 0 2],radius,green)

    ball([0-SHIFT+0.5 0 3],radius,green)
    ball([1+SHIFT+0.5 0 3],radius,green)
    ball([2-SHIFT+0.5 0 3],radius,green)
    ball([3+SHIFT+0.5 0 3],radius,green)
    ball([4-SHIFT+0.5 0 3],radius,green)
    ball([5+SHIFT+0.5 0 3],radius,green)
axis tight
axis equal
view([0 0])
camlight(light)
xlim([-0.5 6])
axis off

disp('wait')
%% A plethora of functions
function ball(c,r,col)
% It just draws a ball
theta=(0:0.01:1)*2*pi;
phi=(0:0.01:1)*pi;
lP=length(phi);
x=r*cos(phi')*sin(theta)+c(1);
y=r*sin(phi')*sin(theta)+c(2);
z=r*ones(lP,1)*cos(theta)+c(3);
surf(x,y,z,'EdgeColor','none','FaceColor',col)
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
