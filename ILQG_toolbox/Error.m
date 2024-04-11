function [error_x,error_y,error_th] = Error(x,xreal,kmax)
error_th = zeros(1,kmax);
error_x = zeros(1,kmax);
error_y = zeros(1,kmax);

for j = 1:kmax
    error_th(j) = norm(x(3,j)-xreal(3,j));
    error_x(j) = norm(x(1,j)-xreal(1,j));
    error_y(j) = norm(x(2,j)-xreal(2,j));
end
error_th = error_th *180/pi;

end
