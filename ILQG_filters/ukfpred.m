function [xest,P] = ukfpred(xest,P,ucorr,dt,Cov_w)

Paug = blkdiag(P,Cov_w); naug = size(Paug,1); dimx = length(xest); dimw = size(Cov_w,1);
alpha_UT = 1; beta_UT = 0; kappa_UT = 3-naug;
[Wm,Wc,ksi,KSI,lambda] = compute_weights(naug,alpha_UT,beta_UT,kappa_UT);

% --- augmented state [theta,x,y,wtheta,wx,wy] ---
xaug = [xest(1);xest(2);xest(3);zeros(dimw,1)];

% ---- generate sigma-points ----
SigPts_01 = [zeros(naug,1) -KSI*eye(naug,naug) KSI*eye(naug,naug)];
SigPts = zeros(naug,2*naug+1);
for j = 1:2*naug+1
    SigPts(:,j) = xaug(:) + sqrtm(Paug)*SigPts_01(:,j);
end
%SigPts = [xaug repmat(xaug,1,naug)-ksi*sqrtm(Paug) repmat(xaug,1,naug)+ksi*sqrtm(Paug)];

% ---- unscented transformation -----
for j = 1:2*naug+1
    xj = SigPts(1:dimx,j);
    wj = SigPts(dimx+1:end,j);
    SigPts(1:dimx,j) = f(xj,ucorr,wj,dt);
end

% ---- mean estimate -----
xest = sum(Wm.*SigPts(1:dimx,:),2);
% ---- covariance -----
P = zeros(dimx,dimx);
for j = 1:2*naug+1
    P = P + Wc(j)*(SigPts(1:dimx,j)-xest)*(SigPts(1:dimx,j)-xest)';
end

end