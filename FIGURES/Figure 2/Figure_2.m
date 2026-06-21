%% Figure 2
% This programme produces Fig. 2 of the paper:
% "Dimensional and Doping Stability of Peierls CDW
%  in Arrays of Coupled One-Dimensional Chains"

%% Clean the worksapce
close all
clear all
clc

%% Set the position of the figure
f=figure;
f.Position=[50,50,[500,1500]];
%% %%%%%%%%%%% %%
%% UPPER PANEL %%
%% %%%%%%%%%%% %%
axes('Position',[ 0.15 0.71 0.8 0.25])
%% Defining constants
N=2;
kappa=0.8;
%% Load data
load(['Electronic_contribution_Half_Filled_N_',...
    num2str(N),'_v2.mat'])
lT=length(tperp);
Epot=4*N*D'.^2*kappa; % Elastic energy
% For different values of tperp,  will have different maps
Nt=1:25:lT;
lNt=length(Nt);
hold off
for T=1:lNt
    c=[0 0 1]*(T-1)/(lNt-1)+[1 0 0]*(lNt-T)/(lNt-1);
    plot(D,(E_electronic(Nt(T),:)+(Epot)'+4/pi*N)*10^3, ...
        'Color',c,'LineWidth',3-1.5*(T-1)/(lNt-1))
    hold on
    plot(-D,(E_electronic(Nt(T),:)+(Epot)'+4/pi*N)*10^3, ...
        'Color',c,'LineWidth',3-1.5*(T-1)/(lNt-1))
end
plot([-0.1 0.1],[0 0],'k--')
ylim([-0.2 0.05])
xlim([-1 1]*0.02)
text(-0.006,-0.14,'Strong coupling','FontSize',13,'Color','b')
text(-0.0058,0.025,'Weak coupling','FontSize',13,'Color','r')
%% Defining constants
t=1; % Hopping parameter
%% Definition of displacement and K/alpha^2
alphaX=(-1:0.00025:1)*0.2*t; % Up to two tenths of the hopping parameter
laX=length(alphaX);
K_over_alpha2=0:0.01:1;
K_over_alpha2_sel=0.5:0.1:1;
lKoA=length(K_over_alpha2);
lKoA_sel=length(K_over_alpha2_sel);
Ee=zeros(lKoA,laX);

yticks(-0.2:0.1:0.1)
box on
xlabel('\Delta (t)')
set(gca,'FontSize',18)


ENERGY=text(-0.0255,-0.17,0,'$E_T- {4t}/{\pi}$ ($t\times 10^{-3}$)','Interpreter','latex','FontSize',15,'Color','k');
set(ENERGY, {'Rotation'}, num2cell(90))
xlabel('$\Delta (t)$','Interpreter','latex')
set(gca,'FontSize',15)
text(-0.027,0.04,0,'(a)','Interpreter','latex','FontSize',20);

%% %%%%%%%%%%%% %%
%% MIDDLE PANEL %%
%% %%%%%%%%%%%% %%
axes('Position',[ 0.15 0.38 0.8 0.27])
N=2;
kappa=[0.00001 0.001:0.001:1.5];
lKA=length(kappa);
%% Load data
load(['Electronic_contribution_Half_Filled_N_',...
    num2str(N),'_v2.mat'])
lT=length(tperp);
Epot=4*N*D'.^2*kappa; % Elastic energy
% For different values of tperp, you have different maps
Minima=zeros(lT,lKA);
for TTT=1:lT
    for K=1:lKA
        Etot=(E_electronic(TTT,:)+(Epot(:,K))'+4/pi*N)*10^3;
        [nMin, idxMin, valMin] = localMinima(Etot);
        if nMin==1 && idxMin>1
            Minima(TTT,K)=2;
        elseif nMin==1 && idxMin==1
            Minima(TTT,K)=0;
        else
            Minima(TTT,K)=1;
        end
    end
end
surf(tperp,kappa,Minima','EdgeColor','none')
hold on
for T=1:lNt
    c=[0 0 1]*(T-1)/(lNt-1)+[1 0 0]*(lNt-T)/(lNt-1);
    if T<3
    scatter3(tperp(Nt(T)),0.8,2,'o', ...
        'MarkerEdgeColor','k', ...
        'MarkerFaceColor',c,'LineWidth',1)
    else
    scatter3(tperp(Nt(T)),0.8,2,'o', ...
        'MarkerEdgeColor','w', ...
        'MarkerFaceColor',c,'LineWidth',1)
    end
end
view([0 90])
text(0.012,-0.25,0,'$\gamma$ (t)','Interpreter','latex','FontSize',16);
ylabel('\kappa (t^{-1})')
%% Analytical solution
tperp_an=0.0001:0.0001:0.025;
lT_an=length(tperp_an);
kappa_an=1/(2*pi)*(log(cos(tperp_an/4)./sin(tperp_an/4))-cos(tperp_an/2));
plot3(tperp_an,kappa_an,ones(1,lT_an)*2,'r--','LineWidth',1)
ylabel('$\kappa$ ($t^{-1}$)','Interpreter','latex')
set(gca,'FontSize',15)
set(gca, 'Layer', 'top')
box on
grid off
text(-0.0045,1.5,0,'(b)','Interpreter','latex','FontSize',20);
ylim([0 1.5])
red=[1 0 0]*0.8;
green=[0 1 0]*0.5;
plot3([0.01 2.8*10^-3],[1 1]*1.01,[1 1]*2,'color',red,'linewidth',2.5)
plot3([2.8*10^-3 0],[1 1]*1.01,[1 1]*2,'color',green,'linewidth',2.5)
plot3([0 0],[1.019 0.981],[1 1]*2,'color',green,'linewidth',2.5)
plot3([0 5.6*10^-3],[1 1]*0.99,[1 1]*2,'color',green,'linewidth',2.5)
plot3([5.6*10^-3 0.01],[1 1]*0.99,[1 1]*2,'color',red,'linewidth',2.5)
plot3([0.01 0.01],[1.019 0.981],[1 1]*2,'color',red,'linewidth',2.5)
text(0.0082,1.08,2, ...
    '$\leftarrow$','Interpreter', ...
    'latex','Color',red,'FontSize',20)
text(0.0057,1.08,2, ...
    '$\leftarrow$','Interpreter', ...
    'latex','Color',red,'FontSize',20)
text(0.0032,1.08,2, ...
    '$\leftarrow$','Interpreter', ...
    'latex','Color',red,'FontSize',20)
text(0.0007,1.08,2, ...
    '$\leftarrow$','Interpreter', ...
    'latex','Color',green,'FontSize',20)

text(0.0007,0.92,2, ...
    '$\rightarrow$','Interpreter', ...
    'latex','Color',green,'FontSize',20)
text(0.0032,0.92,2, ...
    '$\rightarrow$','Interpreter', ...
    'latex','Color',green,'FontSize',20)
text(0.0057,0.92,2, ...
    '$\rightarrow$','Interpreter', ...
    'latex','Color',red,'FontSize',20)
text(0.0082,0.92,2, ...
    '$\rightarrow$','Interpreter', ...
    'latex','Color',red,'FontSize',20)
axes('Position',[ 0.6 0.53 0.3 0.08])
a=1;
R=0.3;
c='g';
ball([0 0 0],R,c)
hold on
ball([1 0 0]*a,R,c)
ball([2 0 0]*a,R,c)
ball([3 0 0]*a,R,c)
ball([4 0 0]*a,R,c)
ball([5 0 0]*a,R,c)
camlight(light)
axis equal
view([0 0])
axis off
axes('Position',[ 0.6 0.55 0.3 0.08])
ball([0 0 0],R,c)
hold on
ball([1 0 0]*a,R,c)
ball([2 0 0]*a,R,c)
ball([3 0 0]*a,R,c)
ball([4 0 0]*a,R,c)
ball([5 0 0]*a,R,c)
camlight(light)
axis equal
view([0 0])
axis off
axes('Position',[ 0.6 0.37 0.3 0.1])
a=1;
x=0.15*a;
R=0.3;
% CDW order state
ball([0+x 0 -1],R,'r')
hold on
ball([1-x 0 -1]*a,R,'r')
ball([2+x 0 -1]*a,R,'r')
ball([3-x 0 -1]*a,R,'r')
ball([4+x 0 -1]*a,R,'r')
ball([5-x 0 -1]*a,R,'r')
camlight(light)
axis equal
view([0 0])
axis off
axes('Position',[ 0.6 0.39 0.3 0.1])
a=1;
x=0.15*a;
R=0.3;
% CDW order state
ball([0+x 0 -1],R,'r')
hold on
ball([1-x 0 -1]*a,R,'r')
ball([2+x 0 -1]*a,R,'r')
ball([3-x 0 -1]*a,R,'r')
ball([4+x 0 -1]*a,R,'r')
ball([5-x 0 -1]*a,R,'r')
camlight(light)
axis equal
view([0 0])
axis off
%% %%%%%%%%%%% %%
%% LOWER PANEL %%
%% %%%%%%%%%%% %%
axes('Position',[ 0.15 0.05 0.8 0.27])
N=2;
kappa=0:0.001:1.5;
lKA=length(kappa);
%% Load data
load('Electronic_contribution_Finite_Doping_N_2_tperp_25e-3t')
lDo=length(doping);
Epot=4*N*D'.^2*kappa; % Elastic energy
% For different values of tperp, you have different maps
Minima=zeros(lDo,lKA);
for Nd=1:lDo
    flag=0;
    for K=1:lKA
        Etot=(E_electronic(Nd,:)+(Epot(:,K))'+4/pi*N)*10^3;
        [nMin, idxMin, valMin] = localMinima(Etot);
        if nMin==1 && idxMin>1
            Minima(Nd,K)=2;
        elseif nMin>1
            flag=1;
            Minima(Nd,K)=1;
        end
        if nMin==1 && flag==1
            Minima(Nd,K)=0;
        end
    end
end

surf(100*doping,kappa,Minima','EdgeColor','none')
hold on
view([0 90])
text(0.44,-0.22,0,'n (\%)','Interpreter','latex','FontSize',16);
ylabel('\kappa (t^{-1})')
set(gca,'FontSize',15)
set(gca, 'Layer', 'top')
grid off
text(-0.17,1.5,0,'(c)','Interpreter','latex','FontSize',20);
box on

%% Mark for the optimal value of the Fermi level
plot3([1 1]*0.4025,[0 1.01],2*[1 1],'r:','linewidth',1)
axes('Position',[ 0.69 0.185 0.23 0.13])
BandStructure_Npar_1Dchains(2,0.025)
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';
ylabel('Energy ', ...
    'color','w','FontSize',15)
ylim([-1 1]*5)
xlabel('Momentum ', ...
    'color','w','FontSize',15)
xticks([])
yticks([])
%% Mark for the optimal value of the Fermi level
plot([-0.5 0.49],[1 1]*tperp*100,'r:','LineWidth',1)
plot([-0.5 0.49],[1 1]*0,'k:','LineWidth',1)
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
