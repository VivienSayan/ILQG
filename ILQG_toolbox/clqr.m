function Lgains = clqr(Q,R,xref,uref,dt)
% u = -Lgains * x

% allocate to store all of the gain matrices
Lgains = ones(size(uref,1),size(xref,1),size(xref,2)); % dim(Lgains) = dim(u)*dim(x)*dim(time)

for t = 1:size(xref,2)
    Aref = [1 0 -dt*uref(1,t)*sin(xref(3,t));...
            0 1 dt*uref(1,t)*cos(xref(3,t));...
            0 0 1                          ];
    Bref = dt*[cos(xref(3,t)) 0;...
               sin(xref(3,t)) 0;...
               0              1];
    [L,S,P] = dlqr(Aref,Bref,Q,R);
    Lgains(:,:,t) = L;
end

end