clear;close all;clc;

TS = 28;            %%% Time steps for neonatal period (days)
TR = 365-TS;        %%% Time steps for the infancy period (days, post infancy)
gS = 1/TS;          %%% grow up from neonates (leave susceptible)
gR = 1/TR;          %%% grow from infancy (leave recovered)

a = 30/1000/TS;     %%% infection rate, S->I, 30 per 1000 sepsis cases
c = 22.34/30/TS;    %%% recovery rate, I->R

d = 22/1000/TS;     %%% death rate not due to sepsis
dI = 7/30/TS;       %%% infected death rate
dR = 48/1000/TR;    %%% non-neonatal death rate for infants
b = 4562;           %%% birth rate, constant input driving susceptibles

T = 1500; %%% number of discrete time steps

%%% Steady state solution
Sinf = b/(d+a+gS);
Iinf = a*b/((d+dI+c)*(d+a+gS));
Rinf = a*b*c/((dR+gR)*(d+dI+c)*(d+a+gS));

%%%% Dynamics, x_k+1 = Ax_k + Bu_k, x = [S;I;R];
A = [1-d-a-gS  0         0;
     a         1-d-dI-c  0;
     0         c         1-dR-gR];
  
%%% Forcing, constant birth rate b
B = [b;0;0];


%%%% INTIAL CONDITION 
x0 = [0;0;0];
x = zeros(3,T);  x(:,1)=x0;
D = zeros(1,T);
for i = 2:T
    
   %%% Discrete time dynamics
   x(:,i) = A*x(:,i-1) + B; %*(1+sin(4*pi*i/365)/6);
   
   %%% Deaths due to sepsis
   D(i) = D(i-1)+dI*x(2,i-1);

end

DT = 1;    %%% discrete time step
t=0:DT:DT*T-DT;



figure;
plot(t,x(1,:)/1000,'color',[0 .5 0],'linewidth',2);hold on;
plot(t,x(2,:)/1000,'r','linewidth',2);hold on;
plot(t,x(3,:)/1000,'b','linewidth',2);hold on;
plot(t,D/1000,'k','linewidth',2);hold on;
plot(t,Sinf*ones(size(t))/1000,'--','color',[0 .5 0],'linewidth',1.5);
plot(t,Iinf*ones(size(t))/1000,'r--','linewidth',1.5);
plot(t,Rinf*ones(size(t))/1000,'b--','linewidth',1.5);
xlabel('Time (days)','fontsize',22);
ylabel('Population (thousands)','fontsize',22);
l = legend('Susceptible, $S_k$','Infected, $I_k$','Recovered, $R_k$','Cumulative Deaths','$S_\infty$','$I_\infty$','$R_\infty$','location','northwest');
set(l,'fontsize',20);
set(l,'interpreter','latex');
set(gca,'fontsize',20);
ylim([0 130]);


