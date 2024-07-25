function [xest,P] = ekfpred(xest,P,ucorr,dt,Cov_w)

% ------ prediction -------
A = [                   1                               0 0;...
     -dt*ucorr(2)*sin(xest(1))-dt*ucorr(3)*cos(xest(1)) 1 0;...
      dt*ucorr(2)*cos(xest(1))-dt*ucorr(3)*sin(xest(1)) 0 1];
B = dt*[1     0            0;...
        0 cos(xest(1)) -sin(xest(1));...
        0 sin(xest(1))  cos(xest(1))];

xest = f(xest,ucorr,[0;0;0],dt); % state propagation

P = A*P*A'+B*Cov_w*B'; % covariance progagation

end