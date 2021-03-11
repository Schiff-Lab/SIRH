clear;close all;clc;

TS = 28;            %%% Time steps for neonatal period (days)
TR = 365-TS;        %%% Time steps for the infancy period (days, post infancy)
gS = 1/TS;          %%% grow up from neonatal class (leave susceptible)
gR = 1/TR;          %%% grow from infancy (leave recovered)

a = 30/1000/TS;     %%% infection rate, S->I, 30 per 1000 sepsis cases
c = 22.34/30/TS;    %%% recovery rate, I->R

d = 22/1000/TS;     %%% death rate not due to sepsis
dI = 7/30/TS;       %%% infected death rate
dR = 47/1000/TR;    %%% non-neonatal death rate for infants not including hydrocephalus
b = 4562;           %%% birth rate, constant input driving susceptibles

h  = 3/22.34/TR;    %%% rate of developing hydrocephalus after recovering
dH = 1/3/TR;        %%% death rate for those that develop hydrocephalus


T = 3000; %%% number of discrete time steps for this simulation

%%% Steady state solution
Sinf = b/(d+a+gS);
Iinf = a*b/((d+dI+c)*(d+a+gS));
Rinf = a*b*c/((dR+gR+h)*(d+dI+c)*(d+a+gS));
Hinf = h*Rinf/(dR+dH);

%%%% Dynamics, x_k+1 = Ax_k + Bu_k, x = [S;I;R];
A = [1-d-a-gS 0         0           0;
     a        1-d-dI-c  0           0;
     0        c         1-dR-gR-h   0;
     0        0         h           1-dR-dH];
  
%%% Forcing, births
B = [b;0;0;0];


%%%% INTIAL CONDITION 
x0 = [0;0;0;0];
x = zeros(4,T);  x(:,1)=x0;
D = zeros(2,T);
for i = 2:T
    
    %%% Discrete time dynamics
    x(:,i) = A*x(:,i-1) + B;%*(1+sin(4*pi*i/365)/6);
    
    %%% Deaths due to infection and hydrocephalus
    D(1,i) = D(1,i-1) + dI*x(2,i-1);
    D(2,i) = D(2,i-1) + dH*x(4,i-1);

end

DT = 1;    %%% discrete time step
t=0:DT:DT*T-DT;



figure;
plot(t,D(1,:)/1000,'r','linewidth',2);hold on;
plot(t,D(2,:)/1000,'k','linewidth',2);hold on;
plot(t,D(1,1501)*ones(size(t))/1000,'r');
plot(t,D(2,2501)*ones(size(t))/1000,'k');
plot(t,D(1,1500+365)*ones(size(t))/1000,'r');
plot(t,D(2,2500+365)*ones(size(t))/1000,'k');
xlabel('Time (days)','fontsize',22);
ylabel('Cumulative Deaths (thousands)','fontsize',22);
l = legend('Sepsis','PIH','1 Year Sepsis','1 Year PIH','location','northwest');
set(l,'fontsize',20);
set(l,'interpreter','latex');
%title('Cumulative Hydrocephalus Deaths','fontsize',24);
set(gca,'fontsize',20);
ylim([0 95]);
text(1500,.9*D(1,1501)/1000,num2str(floor(D(1,1501))),'fontsize',18);
text(1500,1.08*D(1,1500+365)/1000,num2str(floor(D(1,1500+365))),'fontsize',18);
text(2500,.8*D(2,2501)/1000,num2str(floor(D(2,2501))),'fontsize',18);
text(2500,1.2*D(2,2500+365)/1000,num2str(floor(D(2,2500+365))),'fontsize',18);



figure;
plot(t,x(1,:)/1000,'color',[0 .5 0],'linewidth',2);hold on;
plot(t,x(2,:)/1000,'r','linewidth',2);hold on;
plot(t,x(3,:)/1000,'b','linewidth',2);hold on;
plot(t,x(4,:)/1000,'k','linewidth',2);hold on;
plot(t,Sinf*ones(size(t))/1000,'--','color',[0 .5 0],'linewidth',1.5);
plot(t,Iinf*ones(size(t))/1000,'r--','linewidth',1.5);
plot(t,Rinf*ones(size(t))/1000,'b--','linewidth',1.5);
plot(t,Hinf*ones(size(t))/1000,'k--','linewidth',1.5);
xlabel('Time (days)','fontsize',22);
ylabel('Population (thousands)','fontsize',22);
l = legend('Susceptible, $S_k$','Infected, $I_k$','Recovered, $R_k$','Hydrocephalic, $H_k$','location','northwest');
set(l,'fontsize',20);
set(l,'interpreter','latex');
set(gca,'fontsize',20);
ylim([0 130]);


