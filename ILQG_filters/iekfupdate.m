function [xest,P,K] = iekfupdate(xest,P,H,N,z)

% ----- Kalman gain ------
K = P*H'/(H*P*H'+N);
% ----- correction -------
xest = xest + Gamma(xest(3))*K*Rot(-xest(3))*(z-H*xest);
P = P - K*H*P;

end