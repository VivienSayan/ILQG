function [xest,S,P,K] = lukfupdate(xest,S,H,sqrtN,z)

chi = state2chi(xest(3),xest(1:2));

Saug = blkdiag(S,sqrtN); naug = size(Saug,1); dimx = length(xest); dimr = size(sqrtN,1); dimz = length(z);
alpha = 1; beta = 0; kappa = 0;
[Wm,Wc,lambda] = compute_weights(naug,alpha,beta,kappa);
xi = sqrt(naug+lambda);

% --- current augmented state [0;0;0; 0;0;0] ---
xaug = zeros(naug,1);

% ---- generate sigma-points ----
SigPts = xi*[xaug -Saug' Saug'];

% ---- unscented transformation ---
Z = zeros(dimz,2*naug+1);
Z(:,1) = h(chi,zeros(length(sqrtN),1));
for j = 2:2*naug+1
    ksi_j = SigPts(1:dimx,j);
    rj = SigPts(dimx+1:end,j);
    chi_j = chi*expSE2(ksi_j);
    Z(:,j) = h(chi_j,rj);
end

% --- measurement prediction ----
zpred = sum(Wm.*Z,2);
WZ = sqrt(Wc(2:end)).*(Z(:,2:end)-zpred);
[~,Rcz] = qr(WZ');
Sz = Rcz(1:dimz,1:dimz);
Uz = sqrt(abs(Wc(1)))*(Z(:,1) - zpred);
[Sz,~] = cholupdate(Sz,Uz,'-');

% --- cross-covariance ----
Pxz = zeros(dimx,dimz);
for j = 2:2*naug+1
    Pxz  = Pxz + Wc(j)*(SigPts(1:dimx,j))*(Z(:,j)-zpred)';
end

% ---- Kalman gain ----
K = Pxz/Sz/Sz';

% ---- correction ----
ksi_bar = K*(z-zpred);
chi = chi*expSE2(ksi_bar);
[Rot,theta,x] = chi2state(chi);
xest = [x;theta];
U = K*Sz';
for j = 1:dimz
    [S,~] = cholupdate(S,U(:,j),'-');
end 
P = S'*S;

end