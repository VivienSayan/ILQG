function [J,REF_LG,REAL_LG,EST_LG,PEST,UCORR,MEAS] = algo2(filter,xref,uref,xreal0,xest0,P0,M,N,L,Q,R,dt)
% Invariant LQG

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

REF_LG = zeros(3,3,kmax);
xref_lg = zeros(3,kmax);
for t = 1:kmax
    chi = state2chi(xref(3,t),xref(1:2,t));
    REF_LG(:,:,t) = chi;
    [~,theta,x] = chi2state(chi);
    xref_lg(:,t) = [x;theta];
end

REAL_LG = zeros(3,3,kmax); REAL_LG(:,:,1) = state2chi(xreal(3),xreal(1:2));
EST_LG = zeros(3,3,kmax); EST_LG(:,:,1) = state2chi(xest(3),xest(1:2));
PEST = zeros(3,3,kmax); PEST(:,:,1) = P;

Id = invSE2(REF_LG(1:3,1:3,1))*EST_LG(1:3,1:3,1);
dev = logSE2(Id);
%devx = xest(1) - xref(1,1); devy = xest(2) - xref(2,1); devth = atan2(Id(2,1),Id(1,1)); dev = [devx;devy;devth];

switch filter
    case 'iekf'
        ucorr = uref(:,1) - L(:,:,1)*Gamma(-xref(3,1))*(xest-xref(:,1)); 
    case 'left-ukf'
        ucorr = uref(:,1) - L(:,:,1)*dev; 
    case 'left-sr-ukf'
        ucorr = uref(:,1) - L(:,:,1)*dev;
    otherwise
        error('unknown filter')
end

UCORR = zeros(2,kmax); UCORR(:,1) = ucorr; 
MEAS = zeros(2,kmax); MEAS(:,1) = H*xreal + sqrtm(N)*randn(2,1);
KALGAIN = zeros(3,2,kmax); K = KALGAIN(:,:,1);

xbar = xest-xref(:,1);
ubar = ucorr-uref(:,1);
J = xbar'*Q*xbar + ubar'*R*ubar;

for t = 2:kmax

    % ----- real system propagation -----
    m = sqrtm(M)*randn(2,1);
    xreal = f(xreal,ucorr,m,dt);
    REAL_LG(:,:,t) = state2chi(xreal(3),xreal(1:2));

    % ----- prediction --------
    switch filter
        case 'iekf'
            [xest,P] = iekfpred(xest,P,ucorr,dt,M);
        case 'left-ukf' % Left-UKF
            [xest,P] = leftukfpred(xest,P,ucorr,dt,M);
        case'left-sr-ukf' % Left-Square-Root-UKF
            [xest,S,P] = lukfpred(xest,S,ucorr,dt,sqrtM);
        otherwise
            error('unknown filter')
    end

    if obs(t) == 1
        % ----- measurement ------
        z = xreal(1:2) + sqrtm(N)*randn(2,1);
        MEAS(:,t) = z;
    
        % ----- correction -------
        switch filter
            case 'iekf'
                [xest,P,K] = iekfupdate(xest,P,H,N,z);
            case 'left-ukf'
                [xest,P,K] = leftukfupdate(xest,P,H,N,z);
            case 'left-sr-ukf'
                [xest,S,P,K] = lukfupdate(xest,S,H,sqrtN,z);
            otherwise
                error('unknown filter')
        end
    
        EST_LG(:,:,t) = state2chi(xest(3),xest(1:2));
        PEST(:,:,t) = P;
        KALGAIN(:,:,t) = K;
    else
        EST_LG(:,:,t) = state2chi(xest(3),xest(1:2));
        PEST(:,:,t) = P;
        KALGAIN(:,:,t) = K;
    end

    Id = invSE2(REF_LG(1:3,1:3,t))*EST_LG(1:3,1:3,t);
    dev = logSE2(Id);
    %devx = xest(1) - xref(1,t); devy = xest(2) - xref(2,t); devth = atan2(Id(2,1),Id(1,1)); dev = [devx;devy;devth];

    switch filter
        case 'iekf'
            ucorr = uref(:,t) - L(:,:,t)*Gamma(-xref(3,t))*(xest-xref(:,t));
        case 'left-ukf'
            ucorr = uref(:,t) - L(:,:,t)*dev;
        case 'left-sr-ukf'
            ucorr = uref(:,t) - L(:,:,t)*dev;
        otherwise
            error('unknown filter')
    end
    UCORR(:,t) = ucorr; 

    xbar = xest-xref(:,t);
    xbarloc = Gamma(-xest(3))*xbar;
    ubar = ucorr-uref(:,t);
    J = J + xbarloc'*Q*xbarloc + ubar'*R*ubar;
end

end
