function Lt = clqr(Q,R,X,U,dt)
% u(t) = -L(t) * x(t)

% allocation for the gain matrices
Lt = ones(size(U,1),size(X,1),size(X,2)); % dim(Lt) = dim(u)*dim(x)*dim(time)

for t = 1:size(X,2)
    Aref = [                   1                         0 0;...
            -dt*U(2,t)*sin(X(1,t))-dt*U(3,t)*cos(X(1,t)) 1 0;...
            dt*U(2,t)*cos(X(1,t))-dt*U(3,t)*sin(X(1,t))  0 1];
    Bref = dt*[1     0            0;...
               0 cos(X(1,t)) -sin(X(1,t));...
               0 sin(X(1,t))  cos(X(1,t))];
    
    [L,S,e] = dlqr(Aref,Bref,Q,R);
    Lt(:,:,t) = L;
end

end