function [xest,P,K] = ekfupdate(xest,P,H,Cov_v,z)

% ----- Kalman gain ------
K = P*H'/(H*P*H'+Cov_v);
% ----- correction -------
xest = xest + K*(z-H*xest);
P = P - K*H*P;

end