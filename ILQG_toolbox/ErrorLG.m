function [ERROR_TH,ERROR_X,ERROR_Y] = ErrorLG(CHI,CHI_REF,kmax)

ERROR_TH = zeros(1,kmax);
ERROR_X = zeros(1,kmax);
ERROR_Y = zeros(1,kmax);

for j = 1:kmax
    Rot = CHI(1:2,1:2,j);
    Rot_ref = CHI_REF(1:2,1:2,j);
    Id_Error = Rot'*Rot_ref;
    ERROR_TH(j) = norm(atan2(Id_Error(2,1),Id_Error(1,1)));
    ERROR_X(j) = norm(CHI(1,3,j)-CHI_REF(1,3,j));
    ERROR_Y(j) = norm(CHI(2,3,j)-CHI_REF(2,3,j));
end

end