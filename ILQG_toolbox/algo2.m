function [J,XREAL,XEST,PEST,UCORR,MEAS,XREF_LG,XREAL_LG,XEST_LG] = algo2(filter,XREF,UREF,xreal0,xest0,P0,Cov_w_real,Cov_v_real,Cov_w,Cov_v,Lt,Q,R,dt)
% Invariant LQG

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
f_obs = 1; % Hz (measurement frequency)
step = 1/(dt*f_obs);
for i = 1:round(step):t_end
    obs(i) = 1;
end

% convert xref(t) into its equivalent in SE(2)
XREF_LG = zeros(dimx,dimx,t_end);
for t = 1:t_end
    XREF_LG(:,:,t) = state2chi(XREF(1,t),XREF(2:3,t));
end

XREAL = zeros(dimx,t_end); XREAL(:,1) = xreal;
XEST = zeros(dimx,t_end); XEST(:,1) = xest;
XREAL_LG = zeros(dimx,dimx,t_end); XREAL_LG(:,:,1) = state2chi(xreal(1),xreal(2:3));
XEST_LG = zeros(dimx,dimx,t_end); XEST_LG(:,:,1) = state2chi(xest(1),xest(2:3));
PEST = zeros(dimx,dimx,t_end); PEST(:,:,1) = P;
MEAS = zeros(dimz,t_end); MEAS(:,1) = xreal(2:3) + sqrtCov_v_real*randn(dimz,1);
KALMAN_GAIN = zeros(dimx,dimz,t_end); K = KALMAN_GAIN(:,:,1);

switch filter
    case 'iekf'
        ucorr = UREF(:,1) - Lt(:,:,1)*Upsilon(-XREF(1,1)) * (xest-XREF(:,1)); 
    case 'left_ukf'
        Id = invSE2(XREF_LG(1:3,1:3,1)) *XEST_LG(1:3,1:3,1);
        dev = logSE2(Id);
        %devx = xest(2) - XREF(2,1); devy = xest(3) - XREF(3,1); devth = atan2(Id(2,1),Id(1,1)); dev = [devth;devx;devy];
        ucorr = UREF(:,1) - Lt(:,:,1)*dev; 
    case 'srleft_ukf'
        Id = invSE2(XREF_LG(1:3,1:3,1)) *XEST_LG(1:3,1:3,1);
        dev = logSE2(Id);
        %devx = xest(2) - XREF(2,1); devy = xest(3) - XREF(3,1); devth = atan2(Id(2,1),Id(1,1)); dev = [devth;devx;devy];
        ucorr = UREF(:,1) - Lt(:,:,1)*dev;
    otherwise
        error('unknown filter')
end

UCORR = zeros(dimu,t_end); UCORR(:,1) = ucorr; 

xbar = xreal-XREF(:,1); % difference between the true trajectory and the reference trajectory
xbarloc = Upsilon(-xest(1))*xbar;
ubar = ucorr-UREF(:,1); % difference between the corrected input and the reference input
J = xbarloc'*Q*xbarloc + ubar'*R*ubar;

for t = 2:t_end

    % ----- real system propagation -----
    w = sqrtCov_w_real*randn(dimw,1);
    xreal = f(xreal,ucorr,w,dt);
    XREAL(:,t) = xreal;
    XREAL_LG(:,:,t) = state2chi(xreal(1),xreal(2:3));

    % ----- prediction --------
    switch filter
        case 'iekf'
            [xest,P] = iekfpred(xest,P,ucorr,dt,Cov_w);
        case 'left_ukf' % Left-UKF
            [xest,P] = left_ukfpred(xest,P,ucorr,dt,Cov_w);
        case'srleft_ukf' % Square-Root-Left-UKF
            [xest,S,P] = srleft_ukfpred(xest,S,ucorr,dt,sqrtCov_w);
        otherwise
            error('unknown filter')
    end

    if obs(t) == 1
        % ----- measurement ------
        z = xreal(2:3) + sqrtCov_v_real*randn(dimz,1);
        MEAS(:,t) = z;
    
        % ----- correction -------
        switch filter
            case 'iekf'
                [xest,P,K] = iekfupdate(xest,P,H,Cov_v,z);
            case 'left_ukf'
                [xest,P,K] = left_ukfupdate(xest,P,H,Cov_v,z);
            case 'srleft_ukf'
                [xest,S,P,K] = srleft_ukfupdate(xest,S,H,sqrtCov_v,z);
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

    switch filter
        case 'iekf'
            ucorr = UREF(:,t) - Lt(:,:,t)*Upsilon(-XREF(1,t)) * (xest-XREF(:,t));
        case 'left_ukf'
            Id = invSE2(XREF_LG(1:3,1:3,t)) *XEST_LG(1:3,1:3,t);
            dev = logSE2(Id);
            %devx = xest(2) - XREF(2,t); devy = xest(3) - XREF(3,t); devth = atan2(Id(2,1),Id(1,1)); dev = [devth;devx;devy];
            ucorr = UREF(:,t) - Lt(:,:,t)*dev;
        case 'srleft_ukf'
            Id = invSE2(XREF_LG(1:3,1:3,t)) *XEST_LG(1:3,1:3,t);
            [~,theta_err,x_err] = chi2state(Id);
            dev = [theta_err;x_err];
            %dev = logSE2(Id);
            %devx = xest(2) - XREF(2,t); devy = xest(3) - XREF(3,t); devth = atan2(Id(2,1),Id(1,1)); dev = [devth;devx;devy];
            ucorr = UREF(:,t) - Lt(:,:,t)*dev;
        otherwise
            error('unknown filter')
    end
    UCORR(:,t) = ucorr; 

    xbar = xreal-XREF(:,t);
    xbarloc = Upsilon(-xest(1))*xbar;
    ubar = ucorr-UREF(:,t);
    J = J + xbarloc'*Q*xbarloc + ubar'*R*ubar;
end

end
