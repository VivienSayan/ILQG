function [xest,S,P] = QOsrukfpred(xest,S,ucorr,dt,sqrtM)

Saug = blkdiag(S,sqrtM); naug = size(Saug,1); dimx = length(xest); dimq = size(sqrtM,1);
alpha = 1; beta = 0; kappa = 0;
[Wm,Wc,lambda] = compute_weights(naug,alpha,beta,kappa);
xi = sqrt(naug+lambda);

% --- current augmented state [x,y,theta,m1,m2] ---
xaug = [xest(1);xest(2);xest(3);zeros(dimq,1)];

% ---- generate sigma-points ----
SigPts = [xaug repmat(xaug,1,naug)-xi*Saug' repmat(xaug,1,naug)+xi*Saug'];

% ---- optimal quantization -----
mu_th = 1/50; mu_x = 1/10; mu_y = 1/10; mu = diag([mu_th;mu_x;mu_y])*S; N = 300; P = S'*S;
[SigPts,Weights] = QOtheta(mu_th,N,xaug,P,SigPts);

% ---- unscented transformation -----
for j = 1:2*naug+1
    xj = SigPts(1:dimx,j);
    mj = SigPts(dimx+1:end,j);
    SigPts(1:dimx,j) = f(xj,ucorr,mj,dt);
end

% ---- mean estimate -----
xest = sum(Wm.*SigPts(1:dimx,:),2);
% ---- covariance -----
WSigPts = sqrt(Wc(2:end)).*(SigPts(1:dimx,2:end)-xest);
[~,RS] = qr(WSigPts');
S = RS(1:dimx,1:dimx);
Ux = sqrt(abs(Wc(1)))*(SigPts(1:dimx,1) - xest);
[S,~] = cholupdate(S,Ux,'-');
P = S'*S;

end