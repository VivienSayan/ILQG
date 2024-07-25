function xpred = f2(xest,u,w,dt)
% model of the non-holonomic car

% inputs_______________________________
% xpred = [theta; x; y] (rad,m,m)
% u = [vtheta; vx] (rad/s,m/s)
% w = [wtheta; wx] (rad/s, m/s)

% output_________________________
% xpred = [theta; x; y] (rad,m,m)

theta = xest(1) + (u(1)+w(1))*dt;
x     = xest(2) + (u(2)+w(2))*dt* cos(xest(1));
y     = xest(3) + (u(2)+w(2))*dt* sin(xest(1));

xpred = [theta;x;y];

end