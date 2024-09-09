function xpred = f(xest,u,w,dt)
% model of the rigid-body

% inputs_______________________________
% xpred = [theta; x; y] (rad,m,m)
% u = [vtheta; vx; vy] (rad/s,m/s,m/s)
% w = [wtheta; wx; wy] (rad/s, m/s, m/s)

% output_________________________
% xpred = [theta; x; y] (rad,m,m)

eps = 30*pi/180;
if (u(1)+w(1))*dt < eps
    l = (u(2)+w(2))*dt;
    r = (u(3)+w(3))*dt;
else
    l = 1/((u(1)+w(1))*dt) * ((u(2)+w(2))*dt*sin((u(1)+w(1))*dt)-(u(3)+w(3))*dt+(u(3)+w(3))*dt*cos((u(1)+w(1))*dt));
    r = 1/((u(1)+w(1))*dt) * ((u(3)+w(3))*dt*sin((u(1)+w(1))*dt)+(u(2)+w(2))*dt-(u(2)+w(2))*dt*cos((u(1)+w(1))*dt));
end

theta = xest(1) + (u(1)+w(1))*dt;
x     = xest(2) + l*cos(xest(1)) - r*sin(xest(1));
y     = xest(3) + l*sin(xest(1)) + r*cos(xest(1));

xpred = [theta;x;y];

end