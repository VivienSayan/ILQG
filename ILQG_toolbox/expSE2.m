function chi = expSE2(x)

ksi = [x(2);x(3)];
theta = x(1);
if(theta == 0)
    chi = [[eye(2),ksi];...
           [zeros(1,2),1]];
else
    R = [cos(theta) -sin(theta);...
         sin(theta) cos(theta)];
    V = 1/theta*[sin(theta)  , -1+cos(theta);...
                 1-cos(theta), sin(theta)   ];
    Vksi = V*ksi;
    chi = [R Vksi;zeros(1,2) 1];
end

end