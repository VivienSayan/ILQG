function [xest,S,P,K] = sr_ukfupdate(xest,S,H,sqrtCov_v,z)

Saug = blkdiag(S,sqrtCov_v); naug = size(Saug,1); dimx = length(xest); dimv = size(sqrtCov_v,1); dimz = length(z);
alpha_UT = 1; beta_UT = 0; kappa_UT = 3-naug;
[Wm,Wc,ksi,KSI,lambda] = compute_weights(naug,alpha_UT,beta_UT,kappa_UT);

% --- augmented state [theta,x,y,vx,vy] ---
xaug = [xest(1);xest(2);xest(3);zeros(dimv,1)];

% ---- generate sigma-points ----
SigPts_01 = [zeros(naug,1) -KSI*eye(naug,naug) KSI*eye(naug,naug)];
SigPts = zeros(naug,2*naug+1);
for j = 1:2*naug+1
    SigPts(:,j) = xaug(:) + Saug'*SigPts_01(:,j);
end

% --- optimal quantization ---
mu = diag([1/20;1/20;1/20])*S; P = S'*S;
[SigPts,~] = QO(mu(1:1,1:1),100,xaug,P,SigPts,1);

% ---- unscented transformation ---
Z = zeros(dimz,2*naug+1);
for j = 1:2*naug+1
    rj = SigPts(dimx+1:end,j);
    Z(:,j) = SigPts(2:3,j) + rj;
end

% --- measurement prediction ----
zpred = sum(Wm.*Z,2);
WZ = sqrt(Wc(2:end)).*(Z(:,2:end)-zpred);
[~,rSz] = qr(WZ');
Sz = rSz(1:dimz,1:dimz);
Uz = sqrt(abs(Wc(1)))*(Z(:,1) - zpred);
[Sz,~] = cholupdate(Sz,Uz,'-');
Pz = Sz'*Sz;

% --- cross-covariance ----
Pxz = zeros(dimx,dimz);
for j = 1:2*naug+1
    Pxz = Pxz + Wc(j)*(SigPts(1:dimx,j)-xest)*(Z(:,j)-zpred)';
end

% ---- Kalman gain ----
K = Pxz/Pz;

% ---- correction ----
xest = xest + K*(z-zpred);
U = K*Sz';
for j = 1:dimz
    [S,~] = cholupdate(S,U(:,j),'-');
end 
P = S'*S;

end