function [xest,S,P,K] = srukfupdate(xest,S,H,sqrtN,z)

Saug = blkdiag(S,sqrtN); naug = size(Saug,1); dimx = length(xest); dimr = size(sqrtN,1); dimz = length(z);
alpha = 1; beta = 0; kappa = 0;
[Wm,Wc,lambda] = compute_weights(naug,alpha,beta,kappa);
xi = sqrt(naug+lambda);

% --- current augmented state [x,y,theta,n1,n2] ---
xaug = [xest(1);xest(2);xest(3);zeros(dimr,1)];

% ---- generate sigma-points ----
SigPts = [xaug repmat(xaug,1,naug)-xi*Saug' repmat(xaug,1,naug)+xi*Saug'];

% ---- optimal quantization -----
%mu = 1/10*S(1:3,1:3); N = 300; P = S'*S;
%[SigPts,Weights] = QO(mu,N,xaug,P,SigPts);

% ---- unscented transformation ---
Z = zeros(dimz,2*naug+1);
for j = 1:2*naug+1
    rj = SigPts(dimx+1:end,j);
    Z(:,j) = SigPts(1:2,j) + rj;
end

% --- measurement prediction ----
zpred = sum(Wm.*Z,2);
WZ = sqrt(Wc(2:end)).*(Z(:,2:end)-zpred);
[~,RSz] = qr(WZ');
Sz = RSz(1:dimz,1:dimz);
Uz = sqrt(abs(Wc(1)))*(Z(:,1) - zpred);
[Sz,~] = cholupdate(Sz,Uz,'-');

% --- cross-covariance ----
Pxz = zeros(dimx,dimz);
for j = 1:2*naug+1
    Pxz = Pxz + Wc(j)*(SigPts(1:dimx,j)-xest)*(Z(:,j)-zpred)';
end

% ---- Kalman gain ----
Pz = Sz'*Sz;
K = Pxz/Pz;

% ---- correction ----
xest = xest + K*(z-zpred);
U = K*Sz';
for j = 1:dimz
    [S,~] = cholupdate(S,U(:,j),'-');
end 
P = S'*S;

end