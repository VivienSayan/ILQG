function [J,XREAL,XEST,PEST,UCORR,MEAS,XREF_LG,XREAL_LG,XEST_LG] = algo1(filter,XREF,UREF,xreal0,xest0,P0,Cov_w_real,Cov_v_real,Cov_w,Cov_v,Lt,Q,R,dt)
% Conventionnal LQG

t_end = size(XREF,2); dimx = size(XREF,1); dimu = size(UREF,1);
xreal = xreal0;
xest = xest0;
P = P0;
S = chol(P,'upper');

sqrtCov_w_real = chol(Cov_w_real,'upper'); sqrtCov_w = chol(Cov_w,'upper'); dimw = size(Cov_w_real,1);
sqrtCov_v_real = chol(Cov_v_real,'upper'); sqrtCov_v = chol(Cov_v,'upper'); dimv = size(Cov_v_real,1);

% observation matrix y = Hx
H = [0 1 0;...
     0 0 1];
dimz = size(H,1);

obs = zeros(1,t_end);
f_obs = 10; % Hz (measurement frequency)
step = 1/(dt*f_obs);
for i = 1:round(step):t_end
    obs(i) = 1;
end

% convert xref(t) into its equivalent in SE(2)
XREF_LG = zeros(dimx,dimx,t_end);
for t = 1:t_end
    XREF_LG(:,:,t) = state2chi(XREF(1,t),XREF(2:3,t));
end

ucorr = UREF(:,1) - Lt(:,:,1)*(xest-XREF(:,1)); % corrected command input

XREAL = zeros(dimx,t_end); XREAL(:,1) = xreal;
XEST = zeros(dimx,t_end); XEST(:,1) = xest;
XREAL_LG = zeros(dimx,dimx,t_end); XREAL_LG(:,:,1) = state2chi(xreal(1),xreal(2:3));
XEST_LG = zeros(dimx,dimx,t_end); XEST_LG(:,:,1) = state2chi(xest(1),xest(2:3));
PEST = zeros(dimx,dimx,t_end); PEST(:,:,1) = P;
UCORR = zeros(dimu,t_end); UCORR(:,1) = ucorr; 
MEAS = zeros(dimz,t_end); MEAS(:,1) = H*xreal + sqrtCov_v_real*randn(dimz,1);
KALMAN_GAIN = zeros(dimx,dimz,t_end); K = KALMAN_GAIN(:,:,1);

xbar = xest-XREF(:,1); % difference between the true trajectory and the reference trajectory
ubar = ucorr-UREF(:,1); % difference between the corrected input and the reference input
J = xbar'*Q*xbar + ubar'*R*ubar;

for t = 2:t_end

    % ----- real system propagation -----
    w = sqrtCov_w_real*randn(dimw,1);
    xreal = f(xreal,ucorr,w,dt);
    XREAL(:,t) = xreal;
    XREAL_LG(:,:,t) = state2chi(xreal(1),xreal(2:3));

    % ----- prediction --------
    switch filter
        case 'ekf'
            [xest,P] = ekfpred(xest,P,ucorr,dt,Cov_w);
        case 'ukf'
            [xest,P] = ukfpred(xest,P,ucorr,dt,Cov_w);
        case'sr_ukf'
            [xest,S,P] = sr_ukfpred(xest,S,ucorr,dt,sqrtCov_w);
        otherwise
            error('unknown filter')
    end

    if obs(t) == 1
        % ----- measurement ------
        z = xreal(2:3) + sqrtCov_v_real*randn(dimz,1);
        MEAS(:,t) = z;
    
        % ----- correction -------
        switch filter
            case 'ekf'
                [xest,P,K] = ekfupdate(xest,P,H,Cov_v,z);
            case 'ukf'
                [xest,P,K] = ukfupdate(xest,P,H,Cov_v,z);
            case 'sr_ukf'
                [xest,S,P,K] = sr_ukfupdate(xest,S,H,sqrtCov_v,z);
            otherwise
                error('unknown filter')
        end
    
        XEST(:,t) = xest;
        XEST_LG(:,:,t) = state2chi(xest(1),xest(2:3));
        PEST(:,:,t) = P;
        KALMAN_GAIN(:,:,t) = K;
    else
        XEST(:,t) = xest;
        XEST_LG(:,:,t) = state2chi(xest(1),xest(2:3));
        PEST(:,:,t) = P;
        KALMAN_GAIN(:,:,t) = K;
    end

    ucorr = UREF(:,t) - Lt(:,:,t)*(xest-XREF(:,t));
    UCORR(:,t) = ucorr; 

    xbar = xreal-XREF(:,t);
    ubar = ucorr-UREF(:,t);
    J = J + xbar'*Q*xbar + ubar'*R*ubar;
end

end