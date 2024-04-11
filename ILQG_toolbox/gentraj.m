function xref = gentraj(init_state,uref,m,dt)

xref = zeros(3,size(uref,2));
xref(:,1) = init_state;
for t = 2:size(uref,2)
    xref(:,t) = f(xref(:,t-1),uref(:,t-1),m(:,t-1),dt);
end