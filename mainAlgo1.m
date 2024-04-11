%% Invariant Linear Quadratic Gaussian Controller for a simplified car
clear all;
clc; 
close all;

addpath 'ILQG_toolbox';
addpath 'ILQG_filters';
load('ground_truth.mat','u','Tmax','dt','time','kmax','y_GPS');

% set reference command input u (angular velocity and linear velocity)
UREF = zeros(2,kmax);
UREF(1,:) = u(2,:); UREF(1,1:2) = 0.0001;
UREF(2,:) = u(1,:); dimu = size(UREF,1);

% generate (noise free) reference trajectory based on the reference command input
xref0 = [0;0;0]; % initial reference state
m = zeros(2,kmax); % set process noise to zero
XREF = gentraj(xref0,UREF,m,dt); dimx = size(XREF,1);

% compute the L gains via LQR method
Q = [1   0    0;...
      0    1  0;...
      0     0    1]; % state weights

R = [1 0;...
     0 1]; % input weights

% Compute the gain L of the state-feedback u = -Lx that minimizes the the quadratic cost function
Lgains = clqr(Q,R,XREF,UREF,dt);

%%
seed = randi(99999999);
% Kalman parameters
% initial covariance of the estimate
P0 = [(0.3)^2         0          0;...
           0         (0.3)^2     0;...
           0            0      (30*pi/180)^2];
% process noise
M = [(0.1)^2        0;...
        0       (45*pi/180)^2];
% observation noise
N = (0.3)^2*eye(2);


rng(seed);
% real initial state
xreal0 = XREF(:,1) + randn(3,1);
% inital estimate
xest0 = xreal0;%+ sqrtm(P0)*randn(3,1);
[JLQG,XREAL,XEST,PEST,UCORR,MEAS] = algo1('ukf',XREF,UREF,xreal0,xest0,P0,M,N,Lgains,Q,R,dt);

rng(seed);
% real initial state
xreal0 = XREF(:,1) + randn(3,1);
% inital estimate
xest0 = xreal0;%+ sqrtm(P0)*randn(3,1);
[JLQG_QO,XREAL_QO,XEST_QO,PEST_QO,UCORR_QO,MEAS_QO] = algo1('sr_ukf',XREF,UREF,xreal0,xest0,P0,M,N,Lgains,Q,R,dt);

% RESULTS
figure();
subplot(121);
plot(XREF(1,:),XREF(2,:),'k','Linewidth',1); grid on; hold on; xlabel("x(m)"); ylabel("y(m)"); hold on;
plot(MEAS(1,1:1:end),MEAS(2,1:1:end),'.r')
plot(XREAL(1,:),XREAL(2,:),'b')
plot(XEST(1,:),XEST(2,:),'-g')
legend('reference trajectory','measure','real state','estimate')
title('UKF')
subplot(122);
plot(XREF(1,:),XREF(2,:),'k','Linewidth',1); grid on; hold on; xlabel("x(m)"); ylabel("y(m)"); hold on;
plot(MEAS_QO(1,:),MEAS_QO(2,:),'.r')
plot(XREAL_QO(1,:),XREAL_QO(2,:),'b')
plot(XEST_QO(1,:),XEST_QO(2,:),'-g')
legend('reference trajectory','measure','real state','estimate')
title('UKF-Opt')

% TRAJECTORY DEVIATION x-xref
[dev_x,dev_y,dev_th] = Error(XREAL,XREF,kmax); [dev_x_QO,dev_y_QO,dev_th_QO] = Error(XREAL_QO,XREF,kmax);
figure();
subplot(131);
plot(time,dev_x); hold on; plot(time,dev_x_QO);
xlabel('time (s)');
title('|x_{real}-x_{ref}|'); legend('UKF','UKF-Opt');

subplot(132)
plot(time,dev_y); hold on; plot(time,dev_y_QO);
xlabel('time (s)');
title('|y_{real}-y_{ref}|'); legend('UKF','UKF-Opt');

subplot(133)
plot(time,dev_th); hold on; plot(time,dev_th_QO);
xlabel('time (s)');
title('|\theta_{real}-\theta_{ref}|'); legend('UKF','UKF-Opt');

% ESTIMATION ERROR xhat-x
[error_x,error_y,error_th] = Error(XEST,XREAL,kmax); [error_x_QO,error_y_QO,error_th_QO] = Error(XEST_QO,XREAL_QO,kmax);
figure();
subplot(131);
plot(time,error_x); hold on; plot(time,error_x_QO);
xlabel('time (s)');
title('|x_{est}-x_{real}|'); legend('UKF','UKF-Opt');

subplot(132);
plot(time,error_y); hold on; plot(time,error_y_QO);
xlabel('time (s)');
title('|y_{est}-y_{real}|'); legend('UKF','UKF-Opt');

subplot(133);
plot(time,error_th); hold on; plot(time,error_th_QO);
xlabel('time (s)');
title('|\theta_{est}-\theta_{real}|'); legend('UKF','UKF-Opt');

%% STANDARD DEVIATION OF THE ESTIMATE
sigma_x = squeeze(sqrt(PEST(1,1,:)))'; sigma_x_QO = squeeze(sqrt(PEST_QO(1,1,:)))';
sigma_y = squeeze(sqrt(PEST(2,2,:)))'; sigma_y_QO = squeeze(sqrt(PEST_QO(2,2,:)))'; 
sigma_th = squeeze(sqrt(PEST(3,3,:)))'; sigma_th_QO = squeeze(sqrt(PEST_QO(3,3,:)))';

figure();
subplot(131);
plot(time, sigma_x ); hold on; plot(time, sigma_x_QO );
xlabel('t (s)');
ylabel('\sigma_x (m)'); legend('UKF','UKF-Opt');

subplot(132);
plot(time, sigma_y ); hold on; plot(time, sigma_y_QO );
xlabel('t (s)');
ylabel('\sigma_y (m)'); legend('UKF','UKF-Opt');

subplot(133);
plot(time, sigma_th*180/pi ); hold on; plot(time, sigma_th_QO*180/pi );
xlabel('t (s)');
ylabel('\sigma_{\theta} (°)'); legend('UKF','UKF-Opt');

%% CONFIDENCE INTERVAL 3*sigma
figure();
subplot(121)
plot(time,XREAL(1,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST(1,:),'b');
plot(time,XREAL(1,:)+3*sigma_x,'--g');
plot(time,XREAL(1,:)-3*sigma_x,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('x (m)');
subplot(122)
plot(time,XREAL_QO(1,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_QO(1,:),'b');
plot(time,XREAL_QO(1,:)+3*sigma_x_QO,'--g');
plot(time,XREAL_QO(1,:)-3*sigma_x_QO,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('x (m)');

figure();
subplot(121)
plot(time,XREAL(2,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST(2,:),'b');
plot(time,XREAL(2,:)+3*sigma_y,'--g');
plot(time,XREAL(2,:)-3*sigma_y,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('y (m)');
subplot(122)
plot(time,XREAL_QO(2,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_QO(2,:),'b');
plot(time,XREAL_QO(2,:)+3*sigma_y_QO,'--g');
plot(time,XREAL_QO(2,:)-3*sigma_y_QO,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('y (m)');

figure();
subplot(121)
plot(time,XREAL(3,:)*180/pi,'k','linewidth',1); hold on; grid on;
plot(time,XEST(3,:)*180/pi,'b');
plot(time,(XREAL(3,:)+3*sigma_th)*180/pi,'--g');
plot(time,(XREAL(3,:)-3*sigma_th)*180/pi,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('\theta (°)');
subplot(122)
plot(time,XREAL_QO(3,:)*180/pi,'k','linewidth',1); hold on; grid on;
plot(time,XEST_QO(3,:)*180/pi,'b');
plot(time,(XREAL_QO(3,:)+3*sigma_th_QO)*180/pi,'--g');
plot(time,(XREAL_QO(3,:)-3*sigma_th_QO)*180/pi,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('\theta (°)');

%% Animation

figure();

for t = 1:kmax
    plot(xrefLG(1,:),xrefLG(2,:),'k'); grid on; hold on; xlabel("x(m)"); ylabel("y(m)");
    plot(xest(1,:),xest(2,:),'g')
    plot(xreal(1,:),xreal(2,:),'b')

    draw_car(xreal(:,t));
    pause(dt)
    clf;
end

