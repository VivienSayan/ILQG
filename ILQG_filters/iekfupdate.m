function [xest,P,K] = iekfupdate(xest,P,H,Cov_v,z)

% ----- Kalman gain ------
K = P*H'/(H*P*H'+Cov_v);
% ----- correction -------
xest = xest + Upsilon(xest(1))*K*Rot(-xest(1)) *(z-H*xest);
P = P - K*H*P;

end