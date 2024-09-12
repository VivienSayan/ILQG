%% Invariant Linear Quadratic Gaussian Controller for a simplified car
clear all;
clc; 
close all;

addpath 'ILQG_toolbox';
addpath 'ILQG_filters';
load('traj_angle_var_high.mat','u','Tmax','dt','time','kmax','y_GPS','trajReal');
t_end = kmax;

% set reference command input u (angular velocity and linear velocity)
UREF = zeros(3,t_end);
UREF(1,:) = u(1,:); % vtheta (rad/s)
UREF(2,:) = u(2,:); UREF(2,1:2) = 0.0001; % vx (m/s)
UREF(3,:) = u(3,:); % vy (m/s)
dimu = size(UREF,1);

% generate (noise free) reference trajectory based on the reference command input
xref0 = [0;0;0]; % initial reference state [theta;x;y] (in rad,m,m)
W = zeros(dimu,t_end); % set model noise to zero
XREF = gentraj(xref0,UREF,W,dt); dimx = size(XREF,1);

% compute the L gains via LQR method
Q = eye(dimx,dimx); % state weights
R = eye(dimu,dimu); % input weights

% Compute the gain L of the state-feedback u = -Lx that minimizes the the quadratic cost function
Lt = clqr(Q,R,XREF,UREF,dt);

dimz = 2;
%% Simulation
close all; clc;
random_seed = randi(10000);
filter1 = 'ekf';
filter2 = 'sr_ukf';

% real noises
% model noise covariance
Cov_w_real = [(10*pi/180)^2      0        0;...
                0             (0.1)^2     0;...
                0                0     (0.1)^2];
% observation noise covariance
Cov_v_real = (0.5)^2*eye(dimz);

% Kalman parameters-------------------------
% model noise covariance
Cov_w = [(10*pi/180)^2      0         0;...
                0         (0.1)^2     0;...
                0           0      (0.1)^2];
% observation noise covariance
Cov_v = (0.5)^2*eye(dimz);

% initial covariance of the estimate
P0 = [(30*pi/180)^2     0         0;...
           0         (0.3)^2      0;...
           0            0      (0.3)^2];
%-------------------------------------------
Ns = 1;

J_f1_log = zeros(1,Ns);
rng(random_seed);
for n = 1:Ns
% real initial state
xreal0 = XREF(:,1) + randn(3,1);
% inital estimate
xest0 = xreal0 + 0* sqrtm(P0)*randn(3,1);
[JLQG_f1,XREAL_f1,XEST_f1,PEST_f1,UCORR_f1,MEAS_f1,XREF_LG_f1,XREAL_LG_f1,XEST_LG_f1] = algo1(filter1,XREF,UREF,xreal0,xest0,P0,Cov_w_real,Cov_v_real,Cov_w,Cov_v,Lt,Q,R,dt);
J_f1_log(n) = JLQG_f1;
end

J_f2_log = zeros(1,Ns);
rng(random_seed);
for n = 1:Ns
% real initial state
xreal0 = XREF(:,1) + randn(3,1);
% inital estimate
xest0 = xreal0 + 0* sqrtm(P0)*randn(3,1);
[JLQG_f2,XREAL_f2,XEST_f2,PEST_f2,UCORR_f2,MEAS_f2,XREF_LG_f2,XREAL_LG_f2,XEST_LG_f2] = algo1(filter2,XREF,UREF,xreal0,xest0,P0,Cov_w_real,Cov_v_real,Cov_w,Cov_v,Lt,Q,R,dt);
J_f2_log(n) = JLQG_f2;
end

disp('J ekf: ')
mean(J_f1_log)
disp('J ukf: ')
mean(J_f2_log)

std(J_f1_log);
std(J_f2_log);

Jf1f2_mat = [J_f1_log',J_f2_log'];

%% RESULTS
figure();
subplot(121);
plot(XREF(2,:),XREF(3,:),'k','Linewidth',1); grid on; hold on; xlabel("x(m)"); ylabel("y(m)"); hold on;
plot(MEAS_f1(1,1:1:end),MEAS_f1(2,1:1:end),'.r')
plot(XREAL_f1(2,:),XREAL_f1(3,:),'b')
plot(XEST_f1(2,:),XEST_f1(3,:),'-g')
legend('reference trajectory','measure','real state','estimate')
title(filter1)
subplot(122);
plot(XREF(2,:),XREF(3,:),'k','Linewidth',1); grid on; hold on; xlabel("x(m)"); ylabel("y(m)"); hold on;
plot(MEAS_f2(1,:),MEAS_f2(2,:),'.r')
plot(XREAL_f2(2,:),XREAL_f2(3,:),'b')
plot(XEST_f2(2,:),XEST_f2(3,:),'-g')
legend('reference trajectory','measure','real state','estimate')
title(filter1)

% TRAJECTORY DEVIATION xreal-xref
[DEV_TH,DEV_X,DEV_Y] = Error(XREAL_f1,XREF,t_end); [DEV_TH_sr,DEV_X_sr,DEV_Y_sr] = Error(XREAL_f2,XREF,t_end);
figure();
subplot(131);
plot(time,DEV_X); hold on; plot(time,DEV_X_sr);
xlabel('time (s)');
ylabel('|x_{real}-x_{ref}| (m)'); legend(filter1,filter2);

subplot(132)
plot(time,DEV_Y); hold on; plot(time,DEV_Y_sr);
xlabel('time (s)');
ylabel('|y_{real}-y_{ref}| (m)'); legend(filter1,filter2);

subplot(133)
plot(time,DEV_TH*180/pi); hold on; plot(time,DEV_TH_sr*180/pi);
xlabel('time (s)');
ylabel('|\theta_{real}-\theta_{ref}| (°)'); legend(filter1,filter2);

% ESTIMATION ERROR xhat-xreal
[EST_ERROR_TH,EST_ERROR_X,EST_ERROR_Y] = Error(XEST_f1,XREAL_f1,t_end); [EST_ERROR_TH_sr,EST_ERROR_X_sr,EST_ERROR_Y_sr] = Error(XEST_f2,XREAL_f2,t_end);
figure();
subplot(131);
plot(time,EST_ERROR_X); hold on; plot(time,EST_ERROR_X_sr);
xlabel('time (s)');
ylabel('|x_{est}-x_{real}| (m)'); legend(filter1,filter2);

subplot(132);
plot(time,EST_ERROR_Y); hold on; plot(time,EST_ERROR_Y_sr);
xlabel('time (s)');
ylabel('|y_{est}-y_{real}| (m)'); legend(filter1,filter2);

subplot(133);
plot(time,EST_ERROR_TH*180/pi); hold on; plot(time,EST_ERROR_TH_sr*180/pi);
xlabel('time (s)');
ylabel('|\theta_{est}-\theta_{real}| (°)'); legend(filter1,filter2);

%% STANDARD DEVIATION OF THE ESTIMATE
SIGMA_X = squeeze(sqrt(PEST_f1(2,2,:)))'; SIGMA_X_sr = squeeze(sqrt(PEST_f2(2,2,:)))';
SIGMA_Y = squeeze(sqrt(PEST_f1(3,3,:)))'; SIGMA_Y_sr = squeeze(sqrt(PEST_f2(3,3,:)))'; 
SIGMA_TH = squeeze(sqrt(PEST_f1(1,1,:)))'; SIGMA_TH_sr = squeeze(sqrt(PEST_f2(1,1,:)))';

figure();
subplot(131);
plot(time, SIGMA_X ); hold on; 
plot(time, SIGMA_X_sr );
xlabel('t (s)');
ylabel('\sigma_x (m)'); legend(filter1,filter2);

subplot(132);
plot(time, SIGMA_Y ); hold on; 
plot(time, SIGMA_Y_sr );
xlabel('t (s)');
ylabel('\sigma_y (m)'); legend(filter1,filter2);

subplot(133);
plot(time, SIGMA_TH*180/pi ); hold on; 
plot(time, SIGMA_TH_sr*180/pi );
xlabel('t (s)');
ylabel('\sigma_{\theta} (°)'); legend(filter1,filter2);

%% CONFIDENCE INTERVAL 3*sigma
figure();
subplot(121)
plot(time,XREAL_f1(2,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f1(2,:),'b');
plot(time,XREAL_f1(2,:)+3*SIGMA_X,'--g');
plot(time,XREAL_f1(2,:)-3*SIGMA_X,'--g');
title(filter1);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('x (m)');
subplot(122)
plot(time,XREAL_f2(2,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f2(2,:),'b');
plot(time,XREAL_f2(2,:)+3*SIGMA_X_sr,'--g');
plot(time,XREAL_f2(2,:)-3*SIGMA_X_sr,'--g');
title(filter2);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('x (m)');

figure();
subplot(121)
plot(time,XREAL_f1(3,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f1(3,:),'b');
plot(time,XREAL_f1(3,:)+3*SIGMA_Y,'--g');
plot(time,XREAL_f1(3,:)-3*SIGMA_Y,'--g');
title(filter1);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('y (m)');
subplot(122)
plot(time,XREAL_f2(3,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f2(3,:),'b');
plot(time,XREAL_f2(3,:)+3*SIGMA_Y_sr,'--g');
plot(time,XREAL_f2(3,:)-3*SIGMA_Y_sr,'--g');
title(filter2);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('y (m)');

figure();
subplot(121)
plot(time,XREAL_f1(1,:)*180/pi,'k','linewidth',1); hold on; grid on;
plot(time,XEST_f1(1,:)*180/pi,'b');
plot(time,(XREAL_f1(1,:)+3*SIGMA_TH)*180/pi,'--g');
plot(time,(XREAL_f1(1,:)-3*SIGMA_TH)*180/pi,'--g');
title(filter1);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('\theta (°)');
subplot(122)
plot(time,XREAL_f2(1,:)*180/pi,'k','linewidth',1); hold on; grid on;
plot(time,XEST_f2(1,:)*180/pi,'b');
plot(time,(XREAL_f2(1,:)+3*SIGMA_TH_sr)*180/pi,'--g');
plot(time,(XREAL_f2(1,:)-3*SIGMA_TH_sr)*180/pi,'--g');
title(filter2);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('\theta (°)');

%% Animation

figure();
offset_th = 0;
L = 0.1;

for k = 1:t_end
    plot(XREF(2,:),XREF(3,:),'k'); grid on; hold on; xlabel("x(m)"); ylabel("y(m)");
    plot(XEST_f1(2,:),XEST_f1(3,:),'g')
    plot(XREAL_f1(2,:),XREAL_f1(3,:),'b')

    drawRobot(XEST_f1(:,k),offset_th,L);
    pause(dt)
    clf;
end

