function chi = expSE2(xi)

ksi = [xi(1);xi(2)];
theta = xi(3);
if(theta == 0)
    chi = [[eye(2),ksi];...
           [zeros(1,2),1]];
else
    R = [cos(theta) -sin(theta);...
         sin(theta) cos(theta)];
    V = 1/theta*[sin(theta)  , -1+cos(theta);...
                 1-cos(theta), sin(theta)   ];
    Vxi = V*ksi;
    chi = [R Vxi;zeros(1,2) 1];
end

end