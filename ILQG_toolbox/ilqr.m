function Lt = ilqr(Q,R,X,U,dt)
% u(t) = -L(t) * x(t)

dimx = size(Q,1);
dimu = size(R,1);

% allocation for the gain matrices
Lt = ones(size(U,1),size(X,1),size(X,2)); % dim(Lt) = dim(u)*dim(x)*dim(time)

for t = 1:size(X,2)
    Aref = [    1           0            0;...
            -dt*U(3,t)      1         dt*U(1,t);...
             dt*U(2,t)   -dt*U(1,t)      1];
    Bref = dt*eye(dimx,dimu);

    [L,S,e] = dlqr(Aref,Bref,Q,R);
    Lt(:,:,t) = L;
end

end