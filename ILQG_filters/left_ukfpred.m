function [xest,P] = left_ukfpred(xest,P,ucorr,dt,Cov_w)

chi = state2chi(xest(1),xest(2:3));
chi_prev = chi;

Paug = blkdiag(P,Cov_w); naug = size(Paug,1); dimx = length(xest); dimw = size(Cov_w,1);
alpha_UT = 1; beta_UT = 0; kappa_UT = 3-naug;
[Wm,Wc,ksi,KSI,lambda] = compute_weights(naug,alpha_UT,beta_UT,kappa_UT);

% --- current augmented state [0;0;0 ;0;0;0] ---
xaug = zeros(naug,1);

% ---- generate sigma-points ----
SigPts_01 = [zeros(naug,1) -KSI*eye(naug,naug) KSI*eye(naug,naug)];
SigPts = zeros(naug,2*naug+1);
for j = 1:2*naug+1
    SigPts(:,j) = xaug(:) + sqrtm(Paug)*SigPts_01(:,j);
end
%SigPts = ksi*[xaug -sqrtm(Paug) sqrtm(Paug)];

% ----- mean prediction -----
chi = chi*expSE2([ucorr(1);ucorr(2);ucorr(3)]*dt); 
[~,theta,x] = chi2state(chi);
xest = [theta;x];
chi_inv = invSE2(chi);

% ---- unscented transformation -----
for j = 2:2*naug+1
    ksi_j = SigPts(1:dimx,j);
    qj = SigPts(dimx+1:end,j);
    uj = [ucorr(1);ucorr(2);ucorr(3)];
    chi_j = chi_prev*expSE2(ksi_j);
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