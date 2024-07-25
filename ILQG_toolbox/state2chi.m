function chi = state2chi(theta,x)
% embarque theta et x dans matrice SE(2)

Rot = [cos(theta) -sin(theta);...
       sin(theta)  cos(theta)];

chi = [Rot x; zeros(1,2) 1];

end