clear;clc;close all;
% By Donald Ebeigbe and Tyrus Berry

SIRH; %%% define the parameters for the SIRH dynamics

%%% Define the contagious dynamics
beta=0.1/1e5;   %%% Contagious rate
f = @(x) A*x + [-beta*x(1)*x(2);beta*x(1)*x(2);0;0];
Df = @(x) A + [-beta*x(2) -beta*x(1) 0 0;beta*x(2) beta*x(1) 0 0;0 0 0 0;0 0 0 0];

clc;close all; %%% close SIRH plots
rng(2);
N=1500;
W = 2^4*diag([300^2, 25^2, 25^2, 80^2])*1000;  %%% System noise covariance
W = 2^4*diag([300^2, 25^2, 25^2, 80^2])*10;  %%% System noise covariance

c1 = .2/TS;   	%%% Rate that sepsis patients present at care centers
c2 = .6/TR;  	%%% Rate the hydrocephalic patients present at care centers

H = [0 c1 0 0; 0 0 0 c2];   %%% Observation function

x = zeros(4,N+1);           %%% True state trajectory
x(:,1) = (eye(4)-A)\B;      
% for k = 1:1000
%     x(:,1) = f(x(:,1))+B;%%% initial state at equilibrium
% end
Vconstant = diag(H*x(:,1)); %%% Optimal constant choice of V - the true variance at equilibrium
% x(2:4,1) = zeros(3,1);
z = zeros(2,N+1);           %%% Observations
z(:,1) = poissrnd(H*x(:,1));


%%% Initialize the Poisson Kalman Filter (PKF)
xhat = zeros(4,N+1);        %%% State estimate trajectory
xhat(:,1) = x(:,1);         %%% initial state estimate is truth (can add noise)
P = diag([1, 1, 1, 1]);     %%% Covariance estimate

%%% Initialize the standard Kalman Filter (KF)
xhat2 = zeros(4,N+1);        %%% State estimate trajectory
xhat2(:,1) = x(:,1);         %%% initial state estimate is truth (can add noise)
P2 = diag([1, 1, 1, 1]);     %%% Covariance estimate

rtW =sqrtm(W);              %%% Square root of the Q matrix for generating system noise



for k = 1 : N
    
    %B = [b+2000*sin(.05*pi*k/365);0;0;0];  %%% Periodic forcing through birth rate
    
    %%%% Simulate the system, noise, and Poisson measurement %%%%
    
    x(:,k+1) = f(x(:,k)) + B;                                %%% SIRH system
    x(:,k+1) = max(0,x(:,k+1) + rtW*randn(size(x,1),1));    %%% Gaussian system noise
    z(:,k+1) = poissrnd(H*x(:,k+1));                        %%% Poisson observation

    
    
    %%%% Update the PKF %%%%

    % Forecast step
    xhatminus = f(xhat(:,k)) + B;                %%% State forecast
    V = diag(max(.1,H*x(:,k+1)));          %%% V depends on the state estimate, key to the PKF
    V = diag(max(.1,H*xhatminus));
    F = Df(xhat(:,k));
    Pminus = F * P * F' + W;                    %%% Covariance forecast
    K = Pminus * H' / (H * Pminus * H' + V);    %%% Kalman gain

    % Assimilation step
    xhat(:,k+1) = max(0,xhatminus + K*(z(:,k+1) - H*xhatminus));%%% State estimate
    P = (eye(4)-K*H)*Pminus*(eye(4)-K*H)' + K*V*K';             %%% Covariance estimate
    
    
    
    %%%% Update the KF %%%%
    
    V = Vconstant; %%% V is constant as in a standard Kalman filter

    % Forecast step
    xhatminus = f(xhat2(:,k)) + B;               %%% State forecast
    F = Df(xhat2(:,k));
    Pminus = F * P2 * F' + W;                   %%% Covariance forecast
    K = Pminus * H' / (H * Pminus * H' + V);    %%% Kalman gain
    
    % Assimilation step
    xhat2(:,k+1) = max(0,xhatminus + K*(z(:,k+1) - H*xhatminus));%%% State estimate
    P2 = (eye(4)-K*H)*Pminus*(eye(4)-K*H)' + K*V*K';             %%% Covariance estimate
    
    
end


RMSE = sqrt(mean((x - xhat).^2,2))'./std(x')
RMSE2 = sqrt(mean((x - xhat2).^2,2))'./std(x')

big=24;
sml=18;

k = 0 : N;

plotrange = [0 1000];

% Plot the results
close all;
figure;
plot(k, z(1,:), 'or');
hold on;
plot(k, H(1,:)*x, 'k','linewidth',2)
xlabel('Time, k (days)','fontsize',big); 
ylabel('Infected','fontsize',big);
l=legend('Measurement, (y_k)_1', 'Poisson Rate, c_1 I_k');
set(l,'fontsize',sml);
set(gca,'fontsize',sml);
xlim(plotrange);

figure;
plot(k, z(2,:), 'or');
hold on;
plot(k, H(2,:)*x, 'k','linewidth',2)
xlabel('Time, k (days)','fontsize',big); 
ylabel('Hydrocephalic','fontsize',big);
l=legend('Measurement, (y_k)_2', 'Poisson Rate, c_2 H_k');
set(l,'fontsize',sml);
set(gca,'fontsize',sml);
xlim(plotrange);


figure
plot(k, x(1,:), 'k','linewidth',3);hold on;
plot(k, xhat(1,:), '--r','linewidth',1.5);
plot(k, xhat2(1,:), ':b','linewidth',1.5);
plot(k, xhat(1,:), '--r','linewidth',1.5);
ylabel('Susceptible','fontsize',big);
xlabel('Time, k (days)','fontsize',big);
l=legend('Truth, S_k','EPKF','EKF');
set(l,'fontsize',sml);
set(gca,'fontsize',sml);
xlim(plotrange);

figure;
plot(k, x(2,:), 'k','linewidth',3);
hold on;
plot(k, z(1,:)/c1,'go');
plot(k, xhat(2,:), '--r','linewidth',1.5);
plot(k, xhat2(2,:), ':b','linewidth',1.5); 
plot(k, x(2,:), 'k','linewidth',3);
plot(k, xhat(2,:), '--r','linewidth',1.5);
ylabel('Infected','fontsize',big); 
xlabel('Time, k (days)','fontsize',big);
l=legend('Truth, I_k','Rescaled Cases, y_1/c_1','EPKF','EKF');
set(l,'fontsize',sml);
set(gca,'fontsize',sml);
xlim(plotrange);

figure
plot(k, x(3,:), 'k','linewidth',3);hold on;
plot(k, xhat(3,:), '--r','linewidth',1.5);
plot(k, xhat2(3,:), ':b','linewidth',1.5); 
ylabel('Recovered','fontsize',big);
xlabel('Time, k (days)','fontsize',big);
l=legend('Truth, R_k','EPKF','EKF');
set(l,'fontsize',sml);
set(gca,'fontsize',sml);
xlim(plotrange);

figure;
plot(k, x(4,:), 'k','linewidth',3);
hold on;
plot(k, z(2,:)/c2,'go');
plot(k, xhat(4,:), '--r','linewidth',1.5);
plot(k, xhat2(4,:), ':b','linewidth',1.5); 
plot(k, x(4,:), 'k','linewidth',3);
plot(k, xhat(4,:), '--r','linewidth',1.5);
ylabel('Hydrocephalic','fontsize',big); 
xlabel('Time, k (days)','fontsize',big);
l=legend('Truth, H_k','Rescaled Cases, y_2/c_2','EPKF','EKF');
set(l,'fontsize',sml);
set(gca,'fontsize',sml);
xlim(plotrange);



