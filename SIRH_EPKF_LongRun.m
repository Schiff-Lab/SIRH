clear;clc;close all;
% By Donald Ebeigbe and Tyrus Berry

SIRH; %%% define the parameters for the SIRH dynamics

%%% Define the contagious dynamics
beta=0.1/1e5;   %%% Contagious rate
f = @(x) A*x + [-beta*x(1)*x(2);beta*x(1)*x(2);0;0];
Df = @(x) A + [-beta*x(2) -beta*x(1) 0 0;beta*x(2) beta*x(1) 0 0;0 0 0 0;0 0 0 0];

clc;close all; %%% close SIRH plots

Qs = (2.^(-5:0));   %%% multiply the base W by these factors
baseW = 16000*diag([300^2, 25^2, 25^2, 80^2]);

for i=1:length(Qs)

    N=1000000;
    W = Qs(i)*baseW            %%% System noise covariance

    c1 = .2/TS;   	%%% Rate that sepsis patients present at care centers
    c2 = .6/TR;  	%%% Rate the hydrocephalic patients present at care centers

    H = [0 c1 0 0; 0 0 0 c2];   %%% Observation function

    x = zeros(4,N+1);           %%% True state trajectory
    x(:,1) = (eye(4)-A)\B;      
    for k = 1:1000
        x(:,1) = f(x(:,1))+B;%%% initial state at equilibrium
    end
    Vconstant = diag(H*x(:,1)); %%% Optimal constant choice of V - the true variance at equilibrium
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
    
    %%% Initialize the Oracle Poisson Kalman Filter (PKF)
    xhat3 = zeros(4,N+1);        %%% State estimate trajectory
    xhat3(:,1) = x(:,1);         %%% initial state estimate is truth (can add noise)
    P3 = diag([1, 1, 1, 1]);     %%% Covariance estimate

    
    rtW =sqrtm(W);              %%% Square root of the Q matrix for generating system noise


    for k = 1 : N

        %B = [b+2000*sin(.05*pi*k/365);0;0;0];  %%% Periodic forcing through birth rate

        %%%% Simulate the system, noise, and Poisson measurement %%%%

        x(:,k+1) = f(x(:,k)) + B;                                %%% SIRH system
        x(:,k+1) = max(0,x(:,k+1) + rtW*randn(size(x,1),1));    %%% Gaussian system noise
        z(:,k+1) = poissrnd(H*x(:,k+1));                        %%% Poisson observation


        %%%% Update the Oracle PKF %%%%

        % Forecast step
        xhatminus = f(xhat3(:,k)) + B;                %%% State forecast
        V = diag(H*max(.1,x(:,k+1)));          %%% V depends on the state estimate, key to the PKF
        F = Df(xhat3(:,k));
        Pminus = F * P3 * F' + W;                    %%% Covariance forecast
        K = Pminus * H' / (H * Pminus * H' + V);    %%% Kalman gain

        % Assimilation step
        xhat3(:,k+1) = max(0,xhatminus + K*(z(:,k+1) - H*xhatminus));%%% State estimate
        P3 = (eye(4)-K*H)*Pminus*(eye(4)-K*H)' + K*V*K';             %%% Covariance estimate

        
        

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

    RMSE(i,:) = sqrt(mean((x - xhat).^2,2))';
    RMSE2(i,:) = sqrt(mean((x - xhat2).^2,2))';
    RMSE3(i,:) = sqrt(mean((x - xhat3).^2,2))';
    allstd(i,:) = std(x');
    
end



big =22;
sml =18;

figure;
plot(Qs,RMSE3(:,2),'--k','linewidth',2);hold on;
plot(Qs,RMSE(:,2),'r','linewidth',2);
plot(Qs,RMSE2(:,2),'b','linewidth',2);
plot(Qs,RMSE3(:,2),'--k','linewidth',2);
xlabel('System Noise (Multiplier)','fontsize',big); 
ylabel('Infected RMSE','fontsize',big);
l=legend('Oracle EPKF','EPKF','EKF','location','northwest');
set(l,'fontsize',sml);
set(gca,'fontsize',sml);
figure;
plot(Qs,RMSE3(:,4),'--k','linewidth',2);hold on;
plot(Qs,RMSE(:,4),'r','linewidth',2);
plot(Qs,RMSE2(:,4),'b','linewidth',2);
plot(Qs,RMSE3(:,4),'--k','linewidth',2);
xlabel('System Noise (Multiplier)','fontsize',big); 
ylabel('Hydrocephalus RMSE','fontsize',big);
l=legend('Oracle EPKF','EPKF','EKF','location','northwest');
set(l,'fontsize',sml);
set(gca,'fontsize',sml);
    
    
    
