%% Invariant Linear Quadratic Gaussian Controller for a simplified car
clear all;
clc; 
close all;

addpath 'ILQG_toolbox';
addpath 'ILQG_filters';
load('traj_one_loop.mat','u','Tmax','dt','time','kmax','y_GPS','trajReal');
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
Lt = ilqr(Q,R,XREF,UREF,dt);

dimz = 2;
%% Simulation
%clc;
close all;
seed = randi(10000);
filter1 = 'iekf';
filter2 = 'srleft_ukf';
filter3 = 'srleft_ukf_QO';

% real noises-------------------------------
% model noise covariance
Cov_w_real = [(30*pi/180)^2      0         0;...
                0             (0.3)^2      0;...
                0                0        (0.3)^2];
% observation noise covariance
Cov_v_real = (0.2)^2*eye(dimz);

% Kalman parameters-------------------------
% model noise covariance
Cov_w = [(30*pi/180)^2      0          0;...
            0             (0.3)^2      0;...
            0               0        (0.3)^2];
% observation noise covariance
Cov_v = (0.2)^2*eye(dimz);

% ----- initial covariance of the estimate
P0 = [(30*pi/180)^2     0         0;...
           0         (0.3)^2      0;...
           0            0      (0.3)^2];
%-------------------------------------------
Ns = 1;

J_f1_log = zeros(1,Ns);
%rng(seed);
rstream = RandStream('dsfmt19937','Seed',seed);
for n = 1:Ns
    % real initial state
    xreal0 = XREF(:,1) + randn(rstream,3,1);
    % inital estimate
    xest0 = xreal0 + 0* sqrtm(P0)*randn(rstream,3,1);
    [JILQG_f1,XREAL_f1,XEST_f1,PEST_f1,UCORR_f1,MEAS_f1,XREF_LG_f1,XREAL_LG_f1,XEST_LG_f1] = algo2(filter1,XREF,UREF,xreal0,xest0,P0,Cov_w_real,Cov_v_real,Cov_w,Cov_v,Lt,Q,R,dt,rstream);
    J_f1_log(n) = JILQG_f1;
end

J_f2_log = zeros(1,Ns);
%rng(seed);
rstream = RandStream('dsfmt19937','Seed',seed);
for n = 1:Ns
    % real initial state
    xreal0 = XREF(:,1) + randn(rstream,3,1);
    % inital estimate
    xest0 = xreal0+ + 0* sqrtm(P0)*randn(rstream,3,1);
    [JILQG_f2,XREAL_f2,XEST_f2,PEST_f2,UCORR_f2,MEAS_f2,XREF_LG_f2,XREAL_LG_f2,XEST_LG_f2] = algo2(filter2,XREF,UREF,xreal0,xest0,P0,Cov_w_real,Cov_v_real,Cov_w,Cov_v,Lt,Q,R,dt,rstream);
    J_f2_log(n) = JILQG_f2;
end

J_f3_log = zeros(1,Ns);
%rng(seed);
rstream = RandStream('dsfmt19937','Seed',seed);
for n = 1:Ns
    % real initial state
    xreal0 = XREF(:,1) + randn(rstream,3,1);
    % inital estimate
    xest0 = xreal0+ + 0* sqrtm(P0)*randn(rstream,3,1);
    [JILQG_f3,XREAL_f3,XEST_f3,PEST_f3,UCORR_f3,MEAS_f3,XREF_LG_f3,XREAL_LG_f3,XEST_LG_f3] = algo2(filter3,XREF,UREF,xreal0,xest0,P0,Cov_w_real,Cov_v_real,Cov_w,Cov_v,Lt,Q,R,dt,rstream);
    J_f3_log(n) = JILQG_f3;
end
disp('_____')
disp(['J iekf: ',num2str(mean(J_f1_log)),' ± ',num2str(std(J_f1_log))])
disp(['J left-ukf: ',num2str(mean(J_f2_log)),' ± ',num2str(std(J_f2_log))])
disp(['J left-ukf: ',num2str(mean(J_f3_log)),' ± ',num2str(std(J_f3_log))])

%save('tmp.mat');
%% BOXPLOT RESULTS (over the Monte-Carlos)
J_f1_log_clean = rmoutliers(J_f1_log,'percentiles',[0 99]);
J_f2_log_clean = rmoutliers(J_f2_log,'percentiles',[0 99]);
J_f3_log_clean = rmoutliers(J_f3_log,'percentiles',[0 99]);
disp('____')
disp(['J iekf: ',num2str(mean(J_f1_log_clean)),' ± ',num2str(std(J_f1_log_clean))])
disp(['J left-ukf: ',num2str(mean(J_f2_log_clean)),' ± ',num2str(std(J_f2_log_clean))])
disp(['J left-ukf: ',num2str(mean(J_f3_log_clean)),' ± ',num2str(std(J_f3_log_clean))])
Jf1f2f3_mat = [J_f1_log_clean',J_f2_log_clean',J_f3_log_clean'];
figure();
boxplot(Jf1f2f3_mat); hold on; grid on;

%% RESULTS of the last simulation
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
% TRAJECTORY TRACKING
figure();
plot(XREF(2,:),XREF(3,:),'k','Linewidth',1); grid on; hold on; 
plot(MEAS_f1(1,1:1/dt:end),MEAS_f1(2,1:1/dt:end),'.r')
plot(XREAL_f1(2,:),XREAL_f1(3,:),'b','LineWidth',1)
plot(XEST_f1(2,:),XEST_f1(3,:),'--g','LineWidth',1)
legend('reference trajectory','measure','real state','estimate','Interpreter','latex')
xlabel("x(m)",'Interpreter','latex'); ylabel("y(m)",'Interpreter','latex');
title('IEKF','Interpreter','latex')
figure();
plot(XREF(2,:),XREF(3,:),'k','Linewidth',1); grid on; hold on; 
plot(MEAS_f2(1,1:1/dt:end),MEAS_f2(2,1:1/dt:end),'.r')
plot(XREAL_f2(2,:),XREAL_f2(3,:),'b','LineWidth',1)
plot(XEST_f2(2,:),XEST_f2(3,:),'--g','LineWidth',1)
legend('reference trajectory','measure','real state','estimate','Interpreter','latex')
xlabel("x(m)",'Interpreter','latex'); ylabel("y(m)",'Interpreter','latex');
title('Left-UKF-LG','Interpreter','latex')
figure();
plot(XREF(2,:),XREF(3,:),'k','Linewidth',1); grid on; hold on; 
plot(MEAS_f3(1,1:1/dt:end),MEAS_f3(2,1:1/dt:end),'.r')
plot(XREAL_f3(2,:),XREAL_f3(3,:),'b','LineWidth',1)
plot(XEST_f3(2,:),XEST_f3(3,:),'--g','LineWidth',1)
legend('reference trajectory','measure','real state','estimate','Interpreter','latex')
xlabel("x(m)",'Interpreter','latex'); ylabel("y(m)",'Interpreter','latex');
title('OQ-Left-UKF-LG','Interpreter','latex')

% TRAJECTORY TRACKING ERROR xreal-xref
%[DEV_TH_f1,DEV_X_f1,DEV_Y_f1] = Error(XREAL_f1,XREF,t_end); [DEV_TH_f2,DEV_X_f2,DEV_Y_f2] = Error(XREAL_f2,XREF,t_end);
[DEV_TH_f1,DEV_X_f1,DEV_Y_f1] = ErrorLG(XREAL_LG_f1,XREF_LG_f1,t_end); [DEV_TH_f2,DEV_X_f2,DEV_Y_f2] = ErrorLG(XREAL_LG_f2,XREF_LG_f2,t_end); [DEV_TH_f3,DEV_X_f3,DEV_Y_f3] = ErrorLG(XREAL_LG_f3,XREF_LG_f3,t_end);
figure();
subplot(131);
plot(time,DEV_TH_f1*180/pi); hold on;grid on; plot(time,DEV_TH_f2*180/pi); plot(time,DEV_TH_f3*180/pi);
xlabel('time (s)','Interpreter','latex');
ylabel('$|\theta_t-\theta^{*}_t|$ (deg)','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');
subplot(132)
plot(time,DEV_X_f1); hold on;grid on; plot(time,DEV_X_f2); plot(time,DEV_X_f3);
xlabel('time (s)','Interpreter','latex');
ylabel('$|x_t-x^{*}_t|$ (m)','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');
%title('Trajectory tracking error','Interpreter','latex')
subplot(133)
plot(time,DEV_Y_f1); hold on;grid on; plot(time,DEV_Y_f2); plot(time,DEV_Y_f3);
xlabel('time (s)','Interpreter','latex');
ylabel('$|y_t-y^{*}_t|$ (m)','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');


% L2-NORM OF TRAJECTORY TRACKING ERROR
DEV_STATE_f1 = [DEV_TH_f1;DEV_X_f1;DEV_Y_f1]; NORM_DEV_STATE_f1 = vecnorm(DEV_STATE_f1,2,1);
DEV_STATE_f2 = [DEV_TH_f2;DEV_X_f2;DEV_Y_f2]; NORM_DEV_STATE_f2 = vecnorm(DEV_STATE_f2,2,1);
DEV_STATE_f3 = [DEV_TH_f3;DEV_X_f3;DEV_Y_f3]; NORM_DEV_STATE_f3 = vecnorm(DEV_STATE_f3,2,1);
figure()
semilogy(time,NORM_DEV_STATE_f1,'-','LineWidth',1); hold on;grid on; 
semilogy(time,NORM_DEV_STATE_f2,'--','LineWidth',1); 
semilogy(time,NORM_DEV_STATE_f3,'.-','LineWidth',1);
xlabel('time (s)','Interpreter','latex');
ylabel('$||\mathbf{x}_t-\mathbf{x}^{*}_t||_2$','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');
%title('$l_2$-norm of the trajectory tracking error','Interpreter','latex');

% ESTIMATION ERROR xhat-xreal
%[EST_ERROR_TH_f1,EST_ERROR_X_f1,EST_ERROR_Y_f1] = Error(XEST_f1,XREAL_f1,t_end); [EST_ERROR_TH_f2,EST_ERROR_X_f2,EST_ERROR_Y_f2] = Error(XEST_f2,XREAL_f2,t_end);
[EST_ERROR_TH_f1,EST_ERROR_X_f1,EST_ERROR_Y_f1] = ErrorLG(XEST_LG_f1,XREAL_LG_f1,t_end); [EST_ERROR_TH_f2,EST_ERROR_X_f2,EST_ERROR_Y_f2] = ErrorLG(XEST_LG_f2,XREAL_LG_f2,t_end); [EST_ERROR_TH_f3,EST_ERROR_X_f3,EST_ERROR_Y_f3] = ErrorLG(XEST_LG_f3,XREAL_LG_f3,t_end);
figure();
subplot(131);
plot(time,EST_ERROR_TH_f1*180/pi); hold on;grid on; plot(time,EST_ERROR_TH_f2*180/pi); plot(time,EST_ERROR_TH_f3*180/pi);
xlabel('time (s)','Interpreter','latex');
ylabel('$|\hat{\theta}_t-\theta_t|$ (deg)','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');
subplot(132);
plot(time,EST_ERROR_X_f1); hold on;grid on; plot(time,EST_ERROR_X_f2); plot(time,EST_ERROR_X_f3);
xlabel('time (s)','Interpreter','latex');
ylabel('$|\hat{x}_t-x_t|$ (m)','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');
subplot(133);
plot(time,EST_ERROR_Y_f1); hold on;grid on; plot(time,EST_ERROR_Y_f2); plot(time,EST_ERROR_Y_f3);
xlabel('time (s)','Interpreter','latex');
ylabel('$|\hat{y}_t-y_t|$ (m)','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');
%title('Estimation error','Interpreter','latex');

% L2-NORM OF ESTIMATION ERROR
EST_ERROR_STATE_f1 = [EST_ERROR_TH_f1;EST_ERROR_X_f1;DEV_Y_f1]; NORM_EST_ERROR_STATE_f1 = vecnorm(EST_ERROR_STATE_f1,2,1);
EST_ERROR_STATE_f2 = [EST_ERROR_TH_f2;EST_ERROR_X_f2;DEV_Y_f2]; NORM_EST_ERROR_STATE_f2 = vecnorm(EST_ERROR_STATE_f2,2,1);
EST_ERROR_STATE_f3 = [EST_ERROR_TH_f3;EST_ERROR_X_f3;DEV_Y_f3]; NORM_EST_ERROR_STATE_f3 = vecnorm(EST_ERROR_STATE_f3,2,1);
figure()
semilogy(time,NORM_EST_ERROR_STATE_f1,'-','LineWidth',1); hold on;grid on; 
semilogy(time,NORM_EST_ERROR_STATE_f2,'--','LineWidth',1); 
semilogy(time,NORM_EST_ERROR_STATE_f3,'.-','LineWidth',1)
xlabel('time (s)','Interpreter','latex');
ylabel('$||\hat{\mathbf{x}}_t-\mathbf{x}_t||_2$','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');
%title('$l_2$-norm of the state estimation error','Interpreter','latex');


% Norme infinie (erreur d'estimation maximale sur 1 suivi de trajectoire)
norm_inf_f1_est_th = norm(EST_ERROR_TH_f1*180/pi,Inf)
norm_inf_f1_est_x = norm(EST_ERROR_X_f1,Inf)
norm_inf_f1_est_y = norm(EST_ERROR_Y_f1,Inf)

norm_inf_f2_est_th = norm(EST_ERROR_TH_f2*180/pi,Inf)
norm_inf_f2_est_x = norm(EST_ERROR_X_f2,Inf)
norm_inf_f2_est_y = norm(EST_ERROR_Y_f2,Inf)

% Norme infinie (deviation maximale sur 1 suivi de trajectoire)
norm_inf_f1_dev_th = norm(DEV_TH_f1*180/pi,Inf)
norm_inf_f1_dev_x = norm(DEV_X_f1,Inf)
norm_inf_f1_dev_y = norm(DEV_Y_f1,Inf)

norm_inf_f2_dev_th = norm(DEV_TH_f2*180/pi,Inf)
norm_inf_f2_dev_x = norm(DEV_X_f2,Inf)
norm_inf_f2_dev_y = norm(DEV_Y_f2,Inf)

%% STANDARD DEVIATION OF THE ESTIMATE
SIGMA_X = squeeze(sqrt(PEST_f1(2,2,:)))'; SIGMA_X_sr = squeeze(sqrt(PEST_f2(2,2,:)))'; SIGMA_X_srqo = squeeze(sqrt(PEST_f3(2,2,:)))';
SIGMA_Y = squeeze(sqrt(PEST_f1(3,3,:)))'; SIGMA_Y_sr = squeeze(sqrt(PEST_f2(3,3,:)))'; SIGMA_Y_srqo = squeeze(sqrt(PEST_f3(3,3,:)))';
SIGMA_TH = squeeze(sqrt(PEST_f1(1,1,:)))'; SIGMA_TH_sr = squeeze(sqrt(PEST_f2(1,1,:)))'; SIGMA_TH_srqo = squeeze(sqrt(PEST_f3(1,1,:)))';

figure();
subplot(131);
semilogy(time, SIGMA_TH*180/pi ,'-','LineWidth',1); hold on; grid on;
semilogy(time, SIGMA_TH_sr*180/pi ,'--','LineWidth',1);
semilogy(time, SIGMA_TH_srqo*180/pi ,'.-','LineWidth',1);
xlabel('time (sec)','Interpreter','latex');
ylabel('$\sigma_{\theta}$ (deg)','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');
subplot(132);
semilogy(time, SIGMA_X ,'-','LineWidth',1); hold on; grid on;
semilogy(time, SIGMA_X_sr ,'--','LineWidth',1);
semilogy(time, SIGMA_X_srqo ,'.-','LineWidth',1);
xlabel('time (sec)','Interpreter','latex');
ylabel('$\sigma_x$ (m)','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');
%title('Uncertainty of the estimate','Interpreter','latex')
subplot(133);
semilogy(time, SIGMA_Y ,'-','LineWidth',1); hold on; grid on;
semilogy(time, SIGMA_Y_sr ,'--','LineWidth',1);
semilogy(time, SIGMA_Y_srqo ,'.-','LineWidth',1);
xlabel('time (sec)','Interpreter','latex');
ylabel('$\sigma_y$ (m)','Interpreter','latex'); legend('IEKF','Left-UKF-LG','OQ-LEFT-UKF-LG','Interpreter','latex');


%% CRAMER-RAO BOUNDS
Fk_f1 = zeros(dimx,dimx,kmax); Fk_f1(:,:,1) = eye(dimx,dimx);
Fk_f2 = zeros(dimx,dimx,kmax); Fk_f2(:,:,1) = eye(dimx,dimx);
Fk_f3 = zeros(dimx,dimx,kmax); Fk_f3(:,:,1) = eye(dimx,dimx);
for t = 2:kmax
    eps_f1 = logSE2( invSE2( XEST_LG_f1(:,:,t) )*XREAL_LG_f1(:,:,t));
    [Rot_hat, theta_hat, ~] = chi2state(XEST_LG_f1(:,:,t));
    dftrans_deps_f1 = [cos(theta_hat).*dl_deps(eps_f1) - sin(theta_hat).*dr_deps(eps_f1);...
                            sin(theta_hat).*dl_deps(eps_f1) + cos(theta_hat).*dr_deps(eps_f1)];
    Fk_f1(:,:,t) = [1 0 0;...
                   Rot_hat' * dftrans_deps_f1];
%----------
    eps_f2 = logSE2( invSE2( XEST_LG_f2(:,:,t) )*XREAL_LG_f2(:,:,t));
    [Rot_hat, theta_hat, ~] = chi2state(XEST_LG_f2(:,:,t));
    dftrans_deps_f2 = [cos(theta_hat).*dl_deps(eps_f2) - sin(theta_hat).*dr_deps(eps_f2);...
                            sin(theta_hat).*dl_deps(eps_f2) + cos(theta_hat).*dr_deps(eps_f2)];
    Fk_f2(:,:,t) = [1 0 0;...
                    Rot_hat' * dftrans_deps_f2];
%---------
    eps_f3 = logSE2( invSE2( XEST_LG_f3(:,:,t) )*XREAL_LG_f3(:,:,t));
    [Rot_hat, theta_hat, x_hat] = chi2state(XEST_LG_f3(:,:,t));
    dftrans_deps_f3 = [cos(theta_hat).*dl_deps(eps_f3) - sin(theta_hat).*dr_deps(eps_f3);...
                            sin(theta_hat).*dl_deps(eps_f3) + cos(theta_hat).*dr_deps(eps_f3)];
    Fk_f3(:,:,t) = [1 0 0;...
                    Rot_hat' * dftrans_deps_f3];
end

Hk = repmat([0 1 0;0 0 1],1,1,kmax);
Cov_v_real_k = repmat(Cov_v_real,1,1,kmax);
Cov_w_real_k = repmat(Cov_w_real,1,1,kmax);
dimz_k = repmat(dimz,1,kmax);
Bk_f1 = bcr_ilqg(Hk, Fk_f1, Cov_v_real_k, Cov_w_real_k, P0, kmax, dimz_k);
Bk_f2 = bcr_ilqg(Hk, Fk_f2, Cov_v_real_k, Cov_w_real_k, P0, kmax, dimz_k);
Bk_f3 = bcr_ilqg(Hk, Fk_f3, Cov_v_real_k, Cov_w_real_k, P0, kmax, dimz_k);

figure()
subplot(131);
plot(squeeze(Bk_f1(1,1,:)),'k'); hold on; plot(squeeze(Bk_f2(1,1,:)),'--k'); plot(squeeze(Bk_f3(1,1,:)),'.-k');
plot(squeeze(PEST_f1(1,1,:))); plot(squeeze(PEST_f2(1,1,:))); plot(squeeze(PEST_f3(1,1,:)));
xlabel('timestep'); ylabel('P(1,1)')
legend('BCR (IEKF)','BCR (Left-UKF-LG)','BCR (OQ-Left-UKF-LG)','IEKF','Left-UKF-LG','OQ-Left-UKF-LG');

subplot(132);
plot(squeeze(Bk_f1(2,2,:)),'k'); hold on; plot(squeeze(Bk_f2(2,2,:)),'--k'); plot(squeeze(Bk_f3(2,2,:)),'.-k');
plot(squeeze(PEST_f1(2,2,:))); plot(squeeze(PEST_f2(2,2,:))); plot(squeeze(PEST_f3(2,2,:)));
xlabel('timestep'); ylabel('P(2,2)')
legend('BCR (IEKF)','BCR (Left-UKF-LG)','BCR (OQ-Left-UKF-LG)','IEKF','Left-UKF-LG','OQ-Left-UKF-LG');

subplot(133);
plot(squeeze(Bk_f1(3,3,:)),'k'); hold on; plot(squeeze(Bk_f2(3,3,:)),'--k'); plot(squeeze(Bk_f3(3,3,:)),'.-k');
plot(squeeze(PEST_f1(3,3,:))); plot(squeeze(PEST_f2(3,3,:))); plot(squeeze(PEST_f3(3,3,:)));
xlabel('timestep'); ylabel('P(3,3)')
legend('BCR (IEKF)','BCR (Left-UKF-LG)','BCR (OQ-Left-UKF-LG)','IEKF','Left-UKF-LG','OQ-Left-UKF-LG');

%% CONFIDENCE INTERVAL 3*sigma
figure();
subplot(131)
plot(time,XREAL_f1(2,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f1(2,:),'b');
plot(time,XREAL_f1(2,:)+3*SIGMA_X,'--g');
plot(time,XREAL_f1(2,:)-3*SIGMA_X,'--g');
title(filter1);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('x (m)');
subplot(132)
plot(time,XREAL_f2(2,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f2(2,:),'b');
plot(time,XREAL_f2(2,:)+3*SIGMA_X_sr,'--g');
plot(time,XREAL_f2(2,:)-3*SIGMA_X_sr,'--g');
title(filter2);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('x (m)');
subplot(133)
plot(time,XREAL_f3(2,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f3(2,:),'b');
plot(time,XREAL_f3(2,:)+3*SIGMA_X_srqo,'--g');
plot(time,XREAL_f3(2,:)-3*SIGMA_X_srqo,'--g');
title(filter3);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('x (m)');

figure();
subplot(131)
plot(time,XREAL_f1(3,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f1(3,:),'b');
plot(time,XREAL_f1(3,:)+3*SIGMA_Y,'--g');
plot(time,XREAL_f1(3,:)-3*SIGMA_Y,'--g');
title(filter1);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('y (m)');
subplot(132)
plot(time,XREAL_f2(3,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f2(3,:),'b');
plot(time,XREAL_f2(3,:)+3*SIGMA_Y_sr,'--g');
plot(time,XREAL_f2(3,:)-3*SIGMA_Y_sr,'--g');
title(filter2);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('y (m)');
subplot(133)
plot(time,XREAL_f3(3,:),'k','linewidth',1); hold on; grid on;
plot(time,XEST_f3(3,:),'b');
plot(time,XREAL_f3(3,:)+3*SIGMA_Y_srqo,'--g');
plot(time,XREAL_f3(3,:)-3*SIGMA_Y_srqo,'--g');
title(filter3);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('y (m)');

% pour borner de -pi à pi
for tmp = 1:t_end
    [~,XEST_f1(1,tmp),~] = chi2state(XEST_LG_f1(:,:,tmp));
    [~,XREAL_f1(1,tmp),~] = chi2state(XREAL_LG_f1(:,:,tmp));
    [~,XREAL_f2(1,tmp),~] = chi2state(XREAL_LG_f2(:,:,tmp));
    [~,XREAL_f3(1,tmp),~] = chi2state(XREAL_LG_f3(:,:,tmp));
end

figure();
subplot(131)
plot(time,XREAL_f1(1,:)*180/pi,'k','linewidth',1); hold on; grid on;
plot(time,XEST_f1(1,:)*180/pi,'b');
plot(time,(XREAL_f1(1,:)+3*SIGMA_TH)*180/pi,'--g');
plot(time,(XREAL_f1(1,:)-3*SIGMA_TH)*180/pi,'--g');
title(filter1);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('\theta (°)');
subplot(132)
plot(time,XREAL_f2(1,:)*180/pi,'k','linewidth',1); hold on; grid on;
plot(time,XEST_f2(1,:)*180/pi,'b');
plot(time,(XREAL_f2(1,:)+3*SIGMA_TH_sr)*180/pi,'--g');
plot(time,(XREAL_f2(1,:)-3*SIGMA_TH_sr)*180/pi,'--g');
title(filter2);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('\theta (°)');
subplot(133)
plot(time,XREAL_f3(1,:)*180/pi,'k','linewidth',1); hold on; grid on;
plot(time,XEST_f3(1,:)*180/pi,'b');
plot(time,(XREAL_f3(1,:)+3*SIGMA_TH_srqo)*180/pi,'--g');
plot(time,(XREAL_f3(1,:)-3*SIGMA_TH_srqo)*180/pi,'--g');
title(filter3);
legend('real trajectory','estimate','confidence interval 3\sigma');
xlabel('t (s)'); ylabel('\theta (°)');

%% Animation

figure();
offset_th = 0;
L = 0.5;

for k = 1:t_end
    plot(XREF(2,:),XREF(3,:),'k'); grid on; hold on; xlabel("x(m)"); ylabel("y(m)");
    plot(XEST_f2(2,:),XEST_f2(3,:),'g')
    plot(XREAL_f2(2,:),XREAL_f2(3,:),'b')

    drawRobot(XREAL_f2(:,k),offset_th,L);
    pause(dt)
    clf;
end

%%
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
%set(groot, 'DefaultAxesTickLabelInterpreter', 'tex')

figure(879336)
subplot(224)
plot(XREF(2,:),XREF(3,:),'k','Linewidth',0.5); grid on; hold on; 
xlabel("x(m)","Interpreter","latex"); ylabel("y(m)","Interpreter","latex"); hold on;
%legend('reference trajectory','measure','real state','estimate')
title("Scenario 4 - extreme angle variations","Interpreter","latex");
