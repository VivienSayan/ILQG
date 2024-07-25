function [xest,S,P] = sr_ukfpred(xest,S,ucorr,dt,sqrtCov_w)

Saug = blkdiag(S,sqrtCov_w); naug = size(Saug,1); dimx = length(xest); dimw = size(sqrtCov_w,1);
alpha_UT = 1; beta_UT = 0; kappa_UT = 3-naug;
[Wm,Wc,ksi,KSI,lambda] = compute_weights(naug,alpha_UT,beta_UT,kappa_UT);

% --- augmented state [theta,x,y,wtheta,wx,wy] ---
xaug = [xest(1);xest(2);xest(3);zeros(dimw,1)];

% ---- generate sigma-points ----
SigPts_01 = [zeros(naug,1) -KSI*eye(naug,naug) KSI*eye(naug,naug)];
SigPts = zeros(naug,2*naug+1);
for j = 1:2*naug+1
    SigPts(:,j) = xaug(:) + Saug'*SigPts_01(:,j);
end

% --- optimal quantization ---
%mu = diag([1/20;1/20;1/20])*S; P = S'*S;
%[SigPts,~] = QO(mu(1:1,1:1),100,xaug,P,SigPts,1);

% ---- unscented transformation -----
for j = 1:2*naug+1
    xj = SigPts(1:dimx,j);
    wj = SigPts(dimx+1:end,j);
    SigPts(1:dimx,j) = f(xj,ucorr,wj,dt);
end

% ---- mean estimate -----
xest = sum(Wm.*SigPts(1:dimx,:),2);
% ---- covariance -----
WSigPts = sqrt(Wc(2:end)).*(SigPts(1:dimx,2:end)-xest);
[~,rS] = qr(WSigPts');
S = rS(1:dimx,1:dimx);
Ux = sqrt(abs(Wc(1)))*(SigPts(1:dimx,1) - xest);
[S,~] = cholupdate(S,Ux,'-');
P = S'*S;

end