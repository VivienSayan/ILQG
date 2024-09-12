function [xest,P] = iekfpred(xest,P,u,dt,Cov_w)

dimx = length(xest);
dimu = length(u);

w = zeros(dimu,1);
eps = 30*pi/180;
if u(1)*dt < eps
    l = u(2)*dt;
    r = u(3)*dt;
else
    l = 1/((u(1)+w(1))*dt) * ((u(2)+w(2))*dt*sin((u(1)+w(1))*dt)-(u(3)+w(3))*dt+(u(3)+w(3))*dt*cos((u(1)+w(1))*dt));
    r = 1/((u(1)+w(1))*dt) * ((u(3)+w(3))*dt*sin((u(1)+w(1))*dt)+(u(2)+w(2))*dt-(u(2)+w(2))*dt*cos((u(1)+w(1))*dt));
end

% ------ prediction -------
A = [  1       0         0;...
      -r       1      u(1)*dt;...
       l    -u(1)*dt     1];
B = [dt          0  0;...
    -u(3)*dt^2  dt  0;...
     u(2)*dt^2   0 dt];

xest = f(xest,u,zeros(dimu,1),dt); % state propagation

P = A*P*A'+B*Cov_w*B'; % covariance progagation

end