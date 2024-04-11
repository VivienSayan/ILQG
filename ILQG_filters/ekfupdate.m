function [xest,P,K] = ekfupdate(xest,P,H,N,z)

% ----- Kalman gain ------
K = P*H'/(H*P*H'+N);
% ----- correction -------
xest = xest + K*(z-H*xest);
P = P - K*H*P;

end