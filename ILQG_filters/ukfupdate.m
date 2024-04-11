function [xest,P,K] = ukfupdate(xest,P,H,N,z)

Paug = blkdiag(P,N); naug = size(Paug,1); dimx = length(xest); dimr = size(N,1); dimz = length(z);
alpha = 1; beta = 0; kappa = 0;
[Wm,Wc,lambda] = compute_weights(naug,alpha,beta,kappa);
xi = sqrt(naug+lambda);

% --- current augmented state [x,y,theta,n1,n2] ---
xaug = [xest(1);xest(2);xest(3);zeros(dimr,1)];

% ---- generate sigma-points ----
SigPts = [xaug repmat(xaug,1,naug)-xi*sqrtm(Paug) repmat(xaug,1,naug)+xi*sqrtm(Paug)];

% ---- unscented transformation ---
Z = zeros(dimz,2*naug+1);
for j = 1:2*naug+1
    rj = SigPts(dimx+1:end,j);
    Z(:,j) = SigPts(1:2,j) + rj;
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