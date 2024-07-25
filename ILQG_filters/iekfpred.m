function [xest,P] = iekfpred(xest,P,ucorr,dt,Cov_w)

dimx = length(xest);
dimu = length(ucorr);

% ------ prediction -------
A = [   1              0                0;...
     -dt*ucorr(3)      1          dt*ucorr(1);...
      dt*ucorr(2)   -dt*ucorr(1)        1];
B = dt*eye(dimx,dimu);

xest = f(xest,ucorr,zeros(dimu,1),dt); % state propagation

P = A*P*A'+B*Cov_w*B'; % covariance progagation

end