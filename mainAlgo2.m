%% Invariant Linear Quadratic Gaussian Controller for a simplified car
clear all;
clc; 
close all;

addpath 'ILQG_toolbox';
addpath 'ILQG_filters';
addpath 'ILQG_toolboxQO'
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
Linvgains = ilqr(Q,R,XREF,UREF,dt);
%%
seed = randi(999999);
% Kalman parameters
% initial covariance of the estimate
P0 = [(0.35)^2         0          0;...
           0         (0.35)^2     0;...
           0            0      (90*pi/180)^2];
% process noise
M = [(0.01)^2        0;...
        0       (1*pi/180)^2];
% observation noise
N = (0.1)^2*eye(2);

rng(seed);
% real initial state
xreal0 = XREF(:,1) + randn(3,1);
% inital estimate
xest0 = xreal0 + sqrtm(P0)*randn(3,1);
[Jiekf,XREF_G_iekf,XREAL_G_iekf,XEST_G_iekf,PEST_iekf,UCORR_iekf,MEAS_iekf] = algo2('iekf',XREF,UREF,xreal0,xest0,P0,M,N,Linvgains,Q,R,dt);

rng(seed);
% real initial state
xreal0 = XREF(:,1) + randn(3,1);
% inital estimate
xest0 = xreal0 + sqrtm(P0)*randn(3,1);
[Jiukf,XREF_G_iukf,XREAL_G_iukf,XEST_G_iukf,PEST_iukf,UCORR_iukf,MEAS_iukf] = algo2('left-ukf',XREF,UREF,xreal0,xest0,P0,M,N,Linvgains,Q,R,dt);

% Pour tout mettre sur ]-pi;pi]
XREF_iekf = zeros(3,kmax); XREAL_iekf = zeros(3,kmax); XEST_iekf = zeros(3,kmax);
XREF_iukf = zeros(3,kmax); XREAL_iukf = zeros(3,kmax); XEST_iukf = zeros(3,kmax);
for t = 1:kmax
    [~,XREF_iekf(3,t),XREF_iekf(1:2,t)] = chi2state(XREF_G_iekf(:,:,t)); [~,XREF_iukf(3,t),XREF_iukf(1:2,t)] = chi2state(XREF_G_iukf(:,:,t));
    [~,XREAL_iekf(3,t),XREAL_iekf(1:2,t)] = chi2state(XREAL_G_iekf(:,:,t)); [~,XREAL_iukf(3,t),XREAL_iukf(1:2,t)] = chi2state(XREAL_G_iukf(:,:,t));
    [~,XEST_iekf(3,t),XEST_iekf(1:2,t)] = chi2state(XEST_G_iekf(:,:,t)); [~,XEST_iukf(3,t),XEST_iukf(1:2,t)] = chi2state(XEST_G_iukf(:,:,t));
end

% RESULTS
figure();
subplot(121);
plot(XREF_iekf(1,:),XREF_iekf(2,:),'k','Linewidth',1); grid on; hold on; xlabel("x(m)"); ylabel("y(m)"); hold on;
plot(MEAS_iekf(1,:),MEAS_iekf(2,:),'.r')
plot(XREAL_iekf(1,:),XREAL_iekf(2,:),'b')
plot(XEST_iekf(1,:),XEST_iekf(2,:),'g')
legend('reference trajectory','measure','real','estimate')
title('IEKF');
subplot(122);
plot(XREF_iukf(1,:),XREF_iukf(2,:),'k','Linewidth',1); grid on; hold on; xlabel("x(m)"); ylabel("y(m)"); hold on;
plot(MEAS_iukf(1,:),MEAS_iukf(2,:),'.r')
plot(XREAL_iukf(1,:),XREAL_iukf(2,:),'b')
plot(XEST_iukf(1,:),XEST_iukf(2,:),'g')
legend('reference trajectory','measure','real state','estimate')
title('Left-UKF');

% TRAJECTORY DEVIATION x-xref
[dev_x_iekf,dev_y_iekf,dev_th_iekf] = ErrorLG(XREAL_G_iekf,XREF_G_iekf,kmax); [dev_x_iukf,dev_y_iukf,dev_th_iukf] = ErrorLG(XREAL_G_iukf,XREF_G_iukf,kmax);
figure();
subplot(131);
plot(time,dev_x_iekf); hold on; plot(time,dev_x_iukf);
xlabel('time (s)');
title('|x_{real}-x_{ref}|'); legend('IEKF','Left-UKF');

subplot(132)
plot(time,dev_y_iekf); hold on; plot(time,dev_y_iukf);
xlabel('time (s)');
title('|y_{real}-y_{ref}|'); legend('IEKF','Left-UKF');

subplot(133)
plot(time,dev_th_iekf); hold on; plot(time,dev_th_iukf);
xlabel('time (s)');
title('|\theta_{real}-\theta_{ref}|'); legend('IEKF','Left-UKF');

% ESTIMATION ERROR xhat-x
[error_x_iekf,error_y_iekf,error_th_iekf] = ErrorLG(XEST_G_iekf,XREAL_G_iekf,kmax); [error_x_iukf,error_y_iukf,error_th_iukf] = ErrorLG(XEST_G_iukf,XREAL_G_iukf,kmax);
figure();
subplot(131);
plot(time,error_x_iekf); hold on; plot(time,error_x_iukf);
xlabel('time (s)');
title('|x_{est}-x_{real}|'); legend('IEKF','Left-UKF');

subplot(132);
plot(time,error_y_iekf); hold on; plot(time,error_y_iukf);
xlabel('time (s)');
title('|y_{est}-y_{real}|'); legend('IEKF','Left-UKF');

subplot(133);
plot(time,error_th_iekf); hold on; plot(time,error_th_iukf);
xlabel('time (s)');
title('|\theta_{est}-\theta_{real}|'); legend('IEKF','Left-UKF');

%% STANDARD DEVIATION OF THE ESTIMATE
sigma_x_iekf = squeeze(sqrt(PEST_iekf(1,1,:)))'; sigma_x_iukf = squeeze(sqrt(PEST_iukf(1,1,:)))';
sigma_y_iekf = squeeze(sqrt(PEST_iekf(2,2,:)))'; sigma_y_iukf = squeeze(sqrt(PEST_iukf(2,2,:)))'; 
sigma_th_iekf = squeeze(sqrt(PEST_iekf(3,3,:)))'; sigma_th_iukf = squeeze(sqrt(PEST_iukf(3,3,:)))';

figure();
subplot(131);
plot(time, sigma_x_iekf ); hold on; plot(time, sigma_x_iukf );
xlabel('t (s)');
ylabel('\sigma_x (m)'); legend('IEKF','Left-UKF');

subplot(132);
plot(time, sigma_y_iekf ); hold on; plot(time, sigma_y_iukf );
xlabel('t (s)');
ylabel('\sigma_y (m)'); legend('IEKF','Left-UKF');

subplot(133);
plot(time, sigma_th_iekf*180/pi ); hold on; plot(time, sigma_th_iukf*180/pi );
xlabel('t (s)');
ylabel('\sigma_{\theta} (°)'); legend('IEKF','Left-UKF');

%% CONFIDENCE INTERVAL 3*sigma
figure();
subplot(121)
plot(time,XREAL_iekf(1,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_iekf(1,:),'b');
plot(time,XREAL_iekf(1,:)+3*sigma_x_iekf,'--g');
plot(time,XREAL_iekf(1,:)-3*sigma_x_iekf,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('x (m)');
subplot(122)
plot(time,XREAL_iukf(1,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_iukf(1,:),'b');
plot(time,XREAL_iukf(1,:)+3*sigma_x_iukf,'--g');
plot(time,XREAL_iukf(1,:)-3*sigma_x_iukf,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('x (m)');

figure();
subplot(121)
plot(time,XREAL_iekf(2,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_iekf(2,:),'b');
plot(time,XREAL_iekf(2,:)+3*sigma_y_iekf,'--g');
plot(time,XREAL_iekf(2,:)-3*sigma_y_iekf,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('y (m)');
subplot(122)
plot(time,XREAL_iukf(2,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_iukf(2,:),'b');
plot(time,XREAL_iukf(2,:)+3*sigma_y_iukf,'--g');
plot(time,XREAL_iukf(2,:)-3*sigma_y_iukf,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('y (m)');

figure();
subplot(121)
plot(time,XREAL_iekf(3,:)*180/pi,'k','linewidth',1); hold on; grid on;
plot(time,XEST_iekf(3,:)*180/pi,'b');
plot(time,(XREAL_iekf(3,:)+3*sigma_th_iekf)*180/pi,'--g');
plot(time,(XREAL_iekf(3,:)-3*sigma_th_iekf)*180/pi,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('\theta (°)');
subplot(122)
plot(time,XREAL_iukf(3,:)*180/pi,'k','linewidth',1); hold on; grid on;
plot(time,XEST_iukf(3,:)*180/pi,'b');
plot(time,(XEST_iukf(3,:)+3*sigma_th_iukf)*180/pi,'--g');
plot(time,(XREAL_iukf(3,:)-3*sigma_th_iukf)*180/pi,'--g');
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('\theta (°)');

%% Animation

figure();

for t = 1:kmax
    plot(XREF_iekf(1,:),XREF_iekf(2,:),'k'); grid on; hold on; xlabel("x(m)"); ylabel("y(m)");
    plot(XEST_iekf(1,:),XEST_iekf(2,:),'g')
    plot(XREAL_iekf(1,:),XREAL_iekf(2,:),'b')

    draw_car(XREAL_iekf(:,t));
    pause(dt)
    clf;
end

