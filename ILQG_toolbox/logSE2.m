function xi = logSE2(chi)

Rot = chi(1:2,1:2);
theta = atan2(Rot(2,1),Rot(1,1));
x = chi(1:2,3);

if theta == 0
    xi = [x;theta];
else
    A = sin(theta)/theta; B = (1-cos(theta))/theta;
    V_inv = 1/(A^2+B^2)*[A B;-B A];
    
    xi = [V_inv*x;theta];
end

end