%% Figure 1
% This programme produces Fig. 1 of the paper:
% "Dimensional and Doping Stability of Peierls CDW
%  in Arrays of Coupled One-Dimensional Chains"

%% Clean the worksapce
close all
clear all
clc
%% Set the position of the figure
f=figure;
f.Position=[50,50,[700,1500]];
%% %%%%%%%%%%% %%
%% UPPER PANEL %%
%% %%%%%%%%%%% %%
axes('Position',[ 0.15 0.55 0.8 0.4])
%% Defining constants
e=1; % Electron's charge
a=1; % Inter-atomic distance
t=e; % Hopping parameter
%% Defining the momentum mesh
stepK=0.0002;
k=(stepK/2:stepK:(1-stepK/2))*2*pi/(2*a);
dk=k(2)-k(1);
lK=length(k);
%% Definition of displacement and K/alpha^2
alphaX=(-1:0.00025:1)*0.2*t; % Up to two tenths of the hopping parameter
laX=length(alphaX);
K_over_alpha2=0:0.01:1;
K_over_alpha2_sel=0.5:0.1:1;
lKoA=length(K_over_alpha2);
lKoA_sel=length(K_over_alpha2_sel);
Ee=zeros(lKoA,laX);
Ee_an=zeros(lKoA,laX);
for A=1:laX
    %% Analytical approach
    DELTA=alphaX(A)/(2*t);
    Ee_an(:,A)=-4*t/pi*( 1 + ( DELTA )^2/2* ...
                (log( 4/abs(DELTA) ) - 1/2 ));
end
Ee_an_sel=ones(lKoA_sel,1)*Ee_an(1,:);

%% Second term: elastic energy cost
Ep_sel=K_over_alpha2_sel'*(alphaX).^2; % Potential energy of selected rigidities
Etot_sel=Ep_sel+Ee_an_sel; % Total energy of selected rigidities

%% Plot the total energy selected plots
hold on
for KKK=1:lKoA_sel
    c=[0 0 1]*(KKK-1)/(lKoA_sel-1)+[1 0 0]*(lKoA_sel-KKK)/(lKoA_sel-1);
    plot(alphaX/(2*t),(Etot_sel(KKK,:)+4/pi)*10^3,'Color',c,'LineWidth',2)
end
for KKK=1:lKoA_sel
    c=[0 0 1]*(KKK-1)/(lKoA_sel-1)+[1 0 0]*(lKoA_sel-KKK)/(lKoA_sel-1);

    Emin=min(Etot_sel(KKK,:));
    Delta0_sel=4/exp(1)*exp(-2*pi*t*K_over_alpha2_sel(KKK));
    scatter(Delta0_sel,(Emin+4/pi)*10^3,'filled', ...
        'MarkerEdgeColor','k','MarkerFaceColor',c,'LineWidth',2)
    scatter(-Delta0_sel,(Emin+4/pi)*10^3,'filled', ...
        'MarkerEdgeColor','k','MarkerFaceColor',c,'LineWidth',2)
end
view([0 90])
xticks(-0.1:0.05:0.1)
yticks((-1.5:0.5:1.5))
ylim([-1 1]*1.3)
xlim([-0.1 0.1])
plot([-0.1 0.1],[0 0],'k--')
box on
ENERGY=text(-0.12,-0.65,0,'$E_T- 4t/\pi$ ($t\times 10^{-3}$)','Interpreter','latex','FontSize',18,'Color','k');
set(ENERGY, {'Rotation'}, num2cell(90))
xlabel('$\Delta (t)$','Interpreter','latex')
set(gca,'FontSize',15)

text(-0.13,1.2,0,'(a)','Interpreter','latex','FontSize',20);

text(-0.018,0.4,'High rigidity','Color','b','FontSize',15)
text(-0.0175,-0.4,'Low rigidity','Color','r','FontSize',15)

%% %%%%%%%%%%% %%
%% LOWER PANEL %%
%% %%%%%%%%%%% %%
axes('Position',[ 0.15 0.06 0.8 0.4])
%% Defining constants
e=1;
a=1;
t=e; % Hopping parameter
displacement=0.1*a;
K=t/a^2; % Elastic force per unit of atomic displacement
alpha=0.1;
%% Defining the doping levels
stepN=0.0001;
n=(stepN:stepN:0.1)*0.1;
lN=length(n);
%% Definition of displacement and K/alpha^2
alphaX=[0 10^-6]*t; % Up to two tenths of the hopping parameter
laX=length(alphaX);
K_over_alpha2=0.001:0.001:2;
lKoA=length(K_over_alpha2);
Ee_an=zeros(lKoA,laX,lN);
for N=1:lN
    ep=n(N)*pi/2;
    for A=1:laX
        %% Analytical approach
        DELTA=alphaX(A)/(2*t);
        if DELTA==0
            Ee_an(:,A,N)=-4*t/pi*(1-1/2*ep^2);
        else
            Ee_an(:,A,N)=-4*t/pi*(...
                1+DELTA^2/2*(log(4/abs(DELTA))-1/2 )...
                -1/2*(ep*sqrt(ep^2+DELTA^2)+...
                DELTA^2*log((ep+sqrt(ep^2+DELTA^2))/abs(DELTA))) ...
                );
        end
    end
end
%% Second term: elastic energy cost
Ep=K_over_alpha2'*(alphaX).^2;
Etot=Ep+Ee_an;
Diff_E=zeros(lKoA,lN);
for KoA2=1:lKoA
    for N=1:lN
        Diff_E(KoA2,N)=Etot(KoA2,2,N)-Etot(KoA2,1,N);
    end
end
Diff_E=sign(-Diff_E);
%% Plot the phase diagram
surf(n*100,K_over_alpha2,Diff_E,'EdgeColor','none')
view([0 90])
hold on
grid off
doping_CDW=4/(exp(1)*pi)*exp(-2*pi*t*K_over_alpha2);
plot3(doping_CDW*100,K_over_alpha2,ones(1,lKoA),'r--','LineWidth',2)
for KKK=1:lKoA_sel
    c=[0 0 1]*(KKK-1)/(lKoA_sel-1)+[1 0 0]*(lKoA_sel-KKK)/(lKoA_sel-1);
    scatter3(0,K_over_alpha2_sel(KKK),2, ...
        'MarkerEdgeColor','k','MarkerFaceColor',c)
end
xlim([0 0.01*100])
xlabel('$n$ (\%)','Interpreter','latex')
ylabel('$\kappa$ ($t^{-1}$)','Interpreter','latex')
set(gca,'FontSize',15)
set(gca, 'Layer', 'top')
box on
text(-0.14,2,0,'(b)','Interpreter','latex','FontSize',20);
%% Normal state
axes('Position',[ 0.5 0.3 0.4 0.1])
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

axes('Position',[ 0.522 0.33 0.075 0.1])
theta=(0:0.01:1)*pi;
Xc=cos(theta);
Yc=sin(theta);
plot(Xc,Yc,'g','LineWidth',1.5)
axis equal
axis tight
axis off
text(-0.1,1.5,0,'$t$','Interpreter','latex','FontSize',18,'Color','g');
text(-0.3,-1.5,0,'$a$','Interpreter','latex','FontSize',18,'Color','g');

axes('Position',[ 0.5 0.07 0.4 0.1])
a=1;
x=0.15*a;
R=0.3;
%% CDW state
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

axes('Position',[ 0.52 0.097 0.055 0.1])
plot(Xc,Yc,'r','LineWidth',1.5)
axis equal
axis off
axis tight
text(-0.1,1.5,0,'$t_+$','Interpreter','latex','FontSize',18,'Color','r');
text(-1,-2.3,0,'$a-x$','Interpreter','latex','FontSize',13,'Color','r');




axes('Position',[ 0.555 0.098 0.14 0.1])
plot(Xc,Yc,'r','LineWidth',1.5)
axis equal
daspect([1 1.5 1])
axis off
text(-0.1,1.5,0,'$t_-$','Interpreter','latex','FontSize',18,'Color','r');
text(-0.6,-2,0,'$a+x$','Interpreter','latex','FontSize',13,'Color','r');


%% A plethora of functions
function ball(c,r,col)
% It just draws a ball
% theta=(-0.0:0.02:0.5)*2*pi;
% phi=(0:0.05:1)*pi;
theta=(-0.0:0.01:0.5)*2*pi;
phi=(0:0.01:1)*pi;
lP=length(phi);
x=r*cos(phi')*sin(theta)+c(1);
y=r*sin(phi')*sin(theta)+c(2);
z=r*ones(lP,1)*cos(theta)+c(3);
surf(x,y,z,'EdgeColor','none','FaceColor',col)
end