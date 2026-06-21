%% Figure 4
% This programme produces Fig. 4 of the paper:
% "Dimensional and Doping Stability of Peierls CDW
%  in Arrays of Coupled One-Dimensional Chains"

%% Clean the workspace
clear all
close all
clc
%% Set the position of the figure
f=figure;
f.Position=[50,50,[700,1500]];
%% %%%%%%%%%%% %%
%% UPPER PANEL %%
%% %%%%%%%%%%% %%
axes('Position',[ 0.15 0.57 0.8 0.4])
%% Defining constants
kappa=0:0.002:1.5;
lKA=length(kappa);
%% Load data
load('Electronic_contribution_Finite_Doping_skewed_N2_tperp_20e-3t')
Ee_nem=squeeze(Ee_nem);
%% Defining
lD=length(D);
lDo=length(doping);
doping_selected_numbers=43:46:(46*5); 
doping_selected=doping(doping_selected_numbers);
lDo_sel=length(doping_selected);
Epot=4*N*D'.^2*kappa;
% For different values of tperp, you have different maps
Minima=zeros(lDo,lKA);
for Nd=1:lDo
    flag=0;
    for K=1:lKA
        Etot=(Ee_nem(:,Nd)+(Epot(:,K))+4/pi*N)*10^3;
        [nMin, idxMin, valMin] = localMinima(Etot);
        if nMin==1 && idxMin>1
            Minima(Nd,K)=2;
        elseif nMin>1
            flag=0;
            Minima(Nd,K)=2;
        end
        if nMin==1 && flag==1
            Minima(Nd,K)=0;
        end
    end
end
kappa_selected=0.7;
Epots=4*N*D'.^2*kappa_selected;
Etots=zeros(1,lD);

for Nd=1:lDo_sel
        c=[0 0 1]*(Nd-1)/(lDo_sel-1)+[1 0 0]*(lDo_sel-Nd)/(lDo_sel-1);
        Etots=(Ee_nem(:,doping_selected_numbers(Nd))+(Epots)+4/pi*N)*10^3;
        plot(D,Etots,'Color',c,'LineWidth',2)
        hold on
        plot(-D,Etots,'Color',c,'LineWidth',2)
end
xlabel('$\Delta (t)$','Interpreter','latex')
ylabel('$E_{tot}-2t/\pi (10^{-3}t)$','Interpreter','latex')
set(gca,'FontSize',20)
yticks([-0.2 -0.1 0 0.2])


ylim([-0.3 0.04])
xlim([-1 1]*0.02)
text(-0.025,0.03,0,'(a)','Interpreter','latex','FontSize',20);
text(-0.004,-0.27,'Low doping','FontSize',20,'Color','r')
text(-0.0041,0.02,'High doping','FontSize',20,'Color','b')
%% %%%%%%%%%%% %%
%% LOWER PANEL %%
%% %%%%%%%%%%% %%
axes('Position',[ 0.15 0.08 0.8 0.4])
surf(100*doping,kappa,Minima','EdgeColor','none')
hold on
plot3([1 1]*0.455,[0 1]*0.88,[1 1]*2,'r:','LineWidth',1)
plot3([1 1]*0,[0 1]*1.5,[1 1]*2,'r:','LineWidth',1)
for DO=1:lDo_sel
    c=[0 0 1]*(DO-1)/(lDo_sel-1)+[1 0 0]*(lDo_sel-DO)/(lDo_sel-1);
    if DO==3 || DO==4
        scatter3(doping(doping_selected_numbers(DO))*100,0.7,2,'o', ...
            'MarkerEdgeColor','w','SizeData',35, ...
            'MarkerFaceColor',c,'LineWidth',1.5)
    else
        scatter3(doping(doping_selected_numbers(DO))*100,0.7,2,'o', ...
            'MarkerEdgeColor','k','SizeData',35, ...
            'MarkerFaceColor',c,'LineWidth',1.5)
    end
end
view([0 90])
xlabel('n $(\%)$','Interpreter','latex')
ylabel('$\kappa (t^{-1})$','Interpreter','latex')
clim([0 2])
set(gca,'FontSize',20)
axis tight
xlim([0 1])
text(-0.15,1.5,0,'(b)','Interpreter','latex','FontSize',20);

axes('Position',[ 0.73 0.27 0.2 0.2])
dispersions_ParNemZzg_FD(N,tperp)
xlabel('Momentum','Color','w')
ylabel('Energy','Color','w')
set(gca,'FontSize',15)
grid on
plot([-0.5 -0.49],[0 0],'r:','LineWidth',1)
plot([-0.5 -0.49],[1 1]*0.0282843,'r:','LineWidth',1)
ax = gca;
ax.XColor = 'w';
ax.YColor = 'w';
xticks([])
yticks([])
%% %%%%%%%%%%%%%% %%
%% METALLIC STATE %%
%% %%%%%%%%%%%%%% %%
radius=0.3;
axes('Position',[ 0.18 0.33 0.3 0.2])
ball([0 0 0],radius,'g')
hold on
ball([1 0 0],radius,'g')
ball([2 0 0],radius,'g')
ball([3 0 0],radius,'g')
ball([4 0 0],radius,'g')
ball([5 0 0],radius,'g')
ball([0.5 0 1],radius,'g')
ball([1.5 0 1],radius,'g')
ball([2.5 0 1],radius,'g')
ball([3.5 0 1],radius,'g')
ball([4.5 0 1],radius,'g')
ball([5.5 0 1],radius,'g')
axis tight
axis equal
view([0 0])
camlight(light)
xlim([-0.5 6])

axis off
%% %%%%%%%%%%%%% %%
%% ZIG-ZAG STATE %%
%% %%%%%%%%%%%%% %%
SHIFT=0.15;
radius=0.3;
axes('Position',[ 0.18 0.02 0.3 0.2])
ball([0+SHIFT 0 0],radius,'r')
hold on
ball([1-SHIFT 0 0],radius,'r')
ball([2+SHIFT 0 0],radius,'r')
ball([3-SHIFT 0 0],radius,'r')
ball([4+SHIFT 0 0],radius,'r')
ball([5-SHIFT 0 0],radius,'r')

ball([0.5+SHIFT 0 1],radius,'r')
ball([1.5-SHIFT 0 1],radius,'r')
ball([2.5+SHIFT 0 1],radius,'r')
ball([3.5-SHIFT 0 1],radius,'r')
ball([4.5+SHIFT 0 1],radius,'r')
ball([5.5-SHIFT 0 1],radius,'r')
axis tight
axis equal
view([0 0])
camlight(light)
xlim([-0.5 6])
axis off
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
