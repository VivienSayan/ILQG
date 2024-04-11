function [xest,P,K] = left_ukfupdate(xest,P,H,N,z)

chi = state2chi(xest(3),xest(1:2));

Paug = blkdiag(P,N); naug = size(Paug,1); dimx = length(xest); dimr = size(N,1); dimz = length(z);
alpha = 1; beta = 0; kappa = 0;
[Wm,Wc,lambda] = compute_weights(naug,alpha,beta,kappa);
xi = sqrt(naug+lambda);

% --- current augmented state [0;0;0; 0;0;0] ---
xaug = zeros(naug,1);

% ---- generate sigma-points ----
SigPts = xi*[xaug -sqrtm(Paug) sqrtm(Paug)];

% ---- unscented transformation ---
Z = zeros(dimz,2*naug+1);
Z(:,1) = h(chi,zeros(length(N),1));
for j = 2:2*naug+1
    xi_j = SigPts(1:dimx,j);
    rj = SigPts(dimx+1:end,j);
    chi_j = chi*expSE2(xi_j);
    Z(:,j) = h(chi_j,rj);
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
[Rot,theta,x] = chi2state(chi);
xest = [x;theta];
P = P - K*Pz*K';

end