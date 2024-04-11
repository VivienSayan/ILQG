function [J,REAL,EST,PEST,UCORR,MEAS] = algo1(filter,xref,uref,xreal0,xest0,P0,M,N,L,Q,R,dt)
% Conventionnal LQG

kmax = size(xref,2); dimx = size(xref,1); dimu = size(uref,1);
xreal = xreal0;
xest = xest0;

P = P0;
S = chol(P,'upper');
sqrtM = chol(M,'upper');
sqrtN = chol(N,'upper');
H = [1 0 0;...
     0 1 0];

obs = zeros(1,kmax);
DeltaK_obs = 10;
for i = 1:DeltaK_obs:kmax
    obs(i) = 1;
end

ucorr = uref(:,1) - L(:,:,1)*(xest-xref(:,1)); % corrected command input

REAL = zeros(3,kmax); REAL(:,1) = xreal;
EST = zeros(3,kmax); EST(:,1) = xest;
PEST = zeros(3,3,kmax); PEST(:,:,1) = P;
UCORR = zeros(2,kmax); UCORR(:,1) = ucorr; 
MEAS = zeros(2,kmax); MEAS(:,1) = H*xreal + sqrtm(N)*randn(2,1);
KALGAIN = zeros(3,2,kmax); K = KALGAIN(:,:,1);

xbar = xest-xref(:,1);
ubar = ucorr-uref(:,1);
J = xbar'*Q*xbar + ubar'*R*ubar;

for t = 2:kmax

    % ----- real system propagation -----
    m = sqrtM*randn(2,1);
    xreal = f(xreal,ucorr,m,dt);
    REAL(:,t) = xreal;

    % ----- prediction --------
    switch filter
        case 'ekf'
            [xest,P] = ekfpred(xest,P,ucorr,dt,M);
        case 'ukf'
            [xest,P] = ukfpred(xest,P,ucorr,dt,M);
        case'sr_ukf'
            [xest,S,P] = sr_ukfpred(xest,S,ucorr,dt,sqrtM);
        otherwise
            error('unknown filter')
    end

    if obs(t) == 1
        % ----- measurement ------
        z = xreal(1:2) + sqrtm(N)*randn(2,1);
        MEAS(:,t) = z;
    
        % ----- correction -------
        switch filter
            case 'ekf'
                [xest,P,K] = ekfupdate(xest,P,H,N,z);
            case 'ukf'
                [xest,P,K] = ukfupdate(xest,P,H,N,z);
            case 'sr_ukf'
                [xest,S,P,K] = sr_ukfupdate(xest,S,H,sqrtN,z);
            otherwise
                error('unknown filter')
        end
    
        EST(:,t) = xest;
        PEST(:,:,t) = P;
        KALGAIN(:,:,t) = K;
    else
        EST(:,t) = xest;
        PEST(:,:,t) = P;
        KALGAIN(:,:,t) = K;
    end

    ucorr = uref(:,t) - L(:,:,t)*(xest-xref(:,t));
    UCORR(:,t) = ucorr; 

    xbar = xest-xref(:,t);
    ubar = ucorr-uref(:,t);
    J = J + xbar'*Q*xbar + ubar'*R*ubar;
end

end