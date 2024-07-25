function [xest,P,K] = left_ukfupdate(xest,P,H,Cov_v,z)

chi = state2chi(xest(1),xest(2:3));

Paug = blkdiag(P,Cov_v); naug = size(Paug,1); dimx = length(xest); dimv = size(Cov_v,1); dimz = length(z);
alpha_UT = 1; beta_UT = 0; kappa_UT = 3-naug;
[Wm,Wc,ksi,KSI,lambda] = compute_weights(naug,alpha_UT,beta_UT,kappa_UT);

% --- current augmented state [0;0;0; ;0;0] ---
xaug = zeros(naug,1);

% ---- generate sigma-points ----
SigPts_01 = [zeros(naug,1) -KSI*eye(naug,naug) KSI*eye(naug,naug)];
SigPts = zeros(naug,2*naug+1);
for j = 1:2*naug+1
    SigPts(:,j) = xaug(:) + sqrtm(Paug)*SigPts_01(:,j);
end
%SigPts = ksi*[xaug -sqrtm(Paug) sqrtm(Paug)];

% ---- unscented transformation ---
Z = zeros(dimz,2*naug+1);
Z(:,1) = h(chi,zeros(dimv,1));
for j = 2:2*naug+1
    ksi_j = SigPts(1:dimx,j);
    vj = SigPts(dimx+1:end,j);
    chi_j = chi*expSE2(ksi_j);
    Z(:,j) = h(chi_j,vj);
end

% --- measurement prediction ----
zpred = sum(Wm.*Z,2);
Pz = zeros(dimz,dimz);
for j = 1:2*naug+1
    Pz = Pz + Wc(j)*(Z(:,j)-zpred)*(Z(:,j)-zpred)';
end

% --- cross-covariance ----
Pxz = zeros(dimx,dimz);
for j = 2:2*naug+1
    Pxz  = Pxz  + Wc(j)*(SigPts(1:dimx,j))*(Z(:,j)-zpred)';
end

% ---- Kalman gain ----
K = Pxz/Pz;

% ---- correction ----
ksi_bar = K*(z-zpred);
chi = chi*expSE2(ksi_bar);
[~,theta,x] = chi2state(chi);
xest = [theta;x];
P = P - K*Pz*K';
J = JacSE2(ksi_bar,'LEFT');
P = J*P*J';

end