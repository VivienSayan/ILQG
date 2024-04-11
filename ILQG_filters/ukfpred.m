function [xest,P] = ukfpred(xest,P,ucorr,dt,M)

Paug = blkdiag(P,M); naug = size(Paug,1); dimx = length(xest); dimq = size(M,1);
alpha = 1; beta = 0; kappa = 0;
[Wm,Wc,lambda] = compute_weights(naug,alpha,beta,kappa);
xi = sqrt(naug+lambda);

% --- current augmented state [x,y,theta,m1,m2] ---
xaug = [xest(1);xest(2);xest(3);zeros(dimq,1)];

% ---- generate sigma-points ----
SigPts = [xaug repmat(xaug,1,naug)-xi*sqrtm(Paug) repmat(xaug,1,naug)+xi*sqrtm(Paug)];

% ---- optimal quantization -----
%mu = 1/10*S(1:3,1:3); N = 300; P = S'*S;
%[SigPts,Weights] = QO(mu,N,xaug,P,SigPts);

% ---- unscented transformation -----
for j = 1:2*naug+1
    xj = SigPts(1:dimx,j);
    mj = SigPts(dimx+1:end,j);
    SigPts(1:dimx,j) = f(xj,ucorr,mj,dt);
end

% ---- mean estimate -----
xest = sum(Wm.*SigPts(1:dimx,:),2);
% ---- covariance -----
P = zeros(dimx,dimx);
for j = 1:2*naug+1
    P = P + Wc(j)*(SigPts(1:dimx,j)-xest)*(SigPts(1:dimx,j)-xest)';
end

end