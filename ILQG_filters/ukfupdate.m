function [xest,P,K] = ukfupdate(xest,P,H,Cov_v,z)

Paug = blkdiag(P,Cov_v); naug = size(Paug,1); dimx = length(xest); dimv = size(Cov_v,1); dimz = length(z);
alpha_UT = 1; beta_UT = 0; kappa_UT = 3-naug;
[Wm,Wc,ksi,KSI,lambda] = compute_weights(naug,alpha_UT,beta_UT,kappa_UT);

% --- augmented state [theta,x,y,vx,vy] ---
xaug = [xest(1);xest(2);xest(3);zeros(dimv,1)];

% ---- generate sigma-points ----
SigPts_01 = [zeros(naug,1) -KSI*eye(naug,naug) KSI*eye(naug,naug)];
SigPts = zeros(naug,2*naug+1);
for j = 1:2*naug+1
    SigPts(:,j) = xaug(:) + sqrtm(Paug)*SigPts_01(:,j);
end

% ---- unscented transformation ---
Z = zeros(dimz,2*naug+1);
for j = 1:2*naug+1
    rj = SigPts(dimx+1:end,j);
    Z(:,j) = SigPts(2:3,j) + rj;
end

% --- measurement prediction ----
zpred = sum(Wm.*Z,2);
Pz = zeros(dimz,dimz);
for j = 1:2*naug+1
    Pz = Pz + Wc(j)*(Z(:,j)-zpred)*(Z(:,j)-zpred)';
end

% --- cross-covariance ----
Pxz = zeros(dimx,dimz);
for j = 1:2*naug+1
    Pxz = Pxz + Wc(j)*(SigPts(1:dimx,j)-xest)*(Z(:,j)-zpred)';
end

% ---- Kalman gain ----
K = Pxz/Pz;

% ---- correction ----
xest = xest + K*(z-zpred);
P = P - K*Pz*K';

end