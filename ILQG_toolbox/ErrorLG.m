function [error_x,error_y,error_th] = ErrorLG(chi,chi_Real,kmax)

error_th = zeros(1,kmax);
error_x = zeros(1,kmax);
error_y = zeros(1,kmax);

for j = 1:kmax
    Rot = chi(1:2,1:2,j);
    RotReal = chi_Real(1:2,1:2,j);
    Id_Error = Rot'*RotReal;
    error_th(j) = norm(atan2(Id_Error(2,1),Id_Error(1,1)));
    error_x(j) = norm(chi(1,3,j)-chi_Real(1,3,j));
    error_y(j) = norm(chi(2,3,j)-chi_Real(2,3,j));
end
error_th = error_th *180/pi;

end