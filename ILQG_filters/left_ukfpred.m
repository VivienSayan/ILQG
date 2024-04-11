function [xest,P] = left_ukfpred(xest,P,ucorr,dt,M)

chi = state2chi(xest(3),xest(1:2));
chiAnt = chi;

M = diag([M(1,1),0.00001,M(2,2)]);

Paug = blkdiag(P,M); naug = size(Paug,1); dimx = length(xest); dimq = size(M,1);
alpha = 1; beta = 0; kappa = 0;
[Wm,Wc,lambda] = compute_weights(naug,alpha,beta,kappa);
xi = sqrt(naug+lambda);

% --- current augmented state [0;0;0 ;0;0;0] ---
xaug = zeros(naug,1);

% ---- generate sigma-points ----
SigPts = xi*[xaug -sqrtm(Paug) sqrtm(Paug)];

% ----- mean prediction -----
chi = chi*expSE2([ucorr(1);0;ucorr(2)]*dt); 
[Rot,theta,x] = chi2state(chi);
xest = [x;theta];
chi_inv = invSE2(chi);

% ---- unscented transformation -----
for j = 2:2*naug+1
    ksi_j = SigPts(1:dimx,j);
    qj = SigPts(dimx+1:end,j);
    uj = [ucorr(1);0;ucorr(2)];
    chi_j = chiAnt*expSE2(ksi_j);
    chi_j = chi_j*expSE2((uj+qj)*dt);
    Xi_j = chi_inv*chi_j;
    SigPts(1:dimx,j) = logSE2(Xi_j);
end

% ---- mean estimate -----
P = zeros(dimx,dimx);
for j = 1:2*naug+1
    P = P + Wc(j)*SigPts(1:dimx,j)*SigPts(1:dimx,j)';
end

end