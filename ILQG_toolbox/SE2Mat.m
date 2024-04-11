function M = SE2Mat(theta,x)

M = [cos(theta) -sin(theta) x(1);...
     sin(theta) cos(theta)  x(2);...
     0              0        1];

end