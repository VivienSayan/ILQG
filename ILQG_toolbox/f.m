function xpred = f(xest,u,w,dt)
% model of the rigid-body

% inputs_______________________________
% xpred = [theta; x; y] (rad,m,m)
% u = [vtheta; vx; vy] (rad/s,m/s,m/s)
% w = [wtheta; wx; wy] (rad/s, m/s, m/s)

% output_________________________
% xpred = [theta; x; y] (rad,m,m)

theta = xest(1) + (u(1)+w(1))*dt;
x     = xest(2) + (u(2)+w(2))*dt* cos(xest(1));%- (u(3)+w(3))*dt* sin(xest(1));
y     = xest(3) + (u(2)+w(2))*dt* sin(xest(1));%+ (u(3)+w(3))*dt* cos(xest(1));

xpred = [theta;x;y];

end