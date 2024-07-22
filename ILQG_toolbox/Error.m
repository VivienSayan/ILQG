function [ERROR_TH,ERROR_X,ERROR_Y] = Error(X,XREF,t_end)
ERROR_TH = zeros(1,t_end);
ERROR_X = zeros(1,t_end);
ERROR_Y = zeros(1,t_end);

for j = 1:t_end
    ERROR_TH(j) = norm(X(1,j)-XREF(1,j));
    ERROR_X(j)  = norm(X(2,j)-XREF(2,j));
    ERROR_Y(j)  = norm(X(3,j)-XREF(3,j));
end

end
