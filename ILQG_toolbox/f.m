function xpred = f(xest,u,m,dt)
% model of the non-holonomic car
% xpred = [x; y; theta]
% u = [linear_vel ; angular_vel]
% m = [linear_noise; angular_noise]

theta = xest(3) + (u(2)+m(2))*dt;
x =  xest(1) + (u(1)+m(1))*dt*cos(theta);
y = xest(2) + (u(1)+m(1))*dt*sin(theta);

xpred = [x;y;theta];

end