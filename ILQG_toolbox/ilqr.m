function Linvgains = ilqr(Q,R,xref,uref,dt)
% u = -Lgains * x

% allocate to store all of the gain matrices
Linvgains = ones(size(uref,1),size(xref,1),size(xref,2)); % dim(Lgains) = dim(u)*dim(x)*dim(time)

for t = 1:size(xref,2)
    Ainvref = [1             dt*uref(2,t)         0;...
               -dt*uref(2,t)      1         dt*uref(1,t);...
               0                  0            1];
    Binvref = dt*[1 0;...
                  0 0;...
                  0 1];
    [Linv,S,P] = dlqr(Ainvref,Binvref,Q,R);
    Linvgains(:,:,t) = Linv;
end

end