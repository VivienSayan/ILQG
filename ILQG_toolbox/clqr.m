function Lt = clqr(Q,R,X,U,dt)
% u(t) = -L(t) * x(t)

dimx = size(Q,1);
dimu = size(R,1);

% allocation for the gain matrices
Lt = ones(size(U,1),size(X,1),size(X,2)); % dim(Lt) = dim(u)*dim(x)*dim(time)
w = zeros(3,1);

for t = 1:size(X,2)
    u = U(:,t);
    eps = 30*pi/180;
    if u(1)*dt < eps
        l = u(2)*dt;
        r = u(3)*dt;
    else
        l = 1/((u(1)+w(1))*dt) * ((u(2)+w(2))*dt*sin((u(1)+w(1))*dt)-(u(3)+w(3))*dt+(u(3)+w(3))*dt*cos((u(1)+w(1))*dt));
        r = 1/((u(1)+w(1))*dt) * ((u(3)+w(3))*dt*sin((u(1)+w(1))*dt)+(u(2)+w(2))*dt-(u(2)+w(2))*dt*cos((u(1)+w(1))*dt));
    end
    Aref = [                   1          0 0;...
            -l*sin(X(1,t))-r*cos(X(1,t))  1 0;...
             l*cos(X(1,t))-r*sin(X(1,t))  0 1];
    Bref = [                  dt                                   0                0;...
            -cos(X(1,t))*u(3)*dt^2-sin(X(1,t))*u(2)*dt^2    cos(X(1,t))*dt   -sin(X(1,t))*dt;...
            -sin(X(1,t))*u(3)*dt^2+cos(X(1,t))*u(2)*dt^2    sin(X(1,t))*dt    cos(X(1,t))*dt];
    
    [L,S,e] = dlqr(Aref,Bref,Q,R);
    Lt(:,:,t) = L;
end

end