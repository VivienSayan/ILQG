function J = JacSE2(ksi,filter) 
% RIGHT jacobian is for filter = LEFT-ukf
% LEFT jacobian is for filter = RIGHT-ukf

th = ksi(1);
v1 = ksi(2);
v2 = ksi(3);

tolerance = 0.00001*pi/180;

if abs(th) < tolerance % If |theta| is small, fall back on the series representation
    switch filter
        case 'LEFT' % <=> RIGHT jacobian
            %J = eye(3);
            J = Jacright_Series(ksi,10);
        case 'RIGHT' % <=> LEFT jacobian
            %J = eye(3);
            J = Jacleft_Series(ksi,10);
        otherwise
            error('filter does not exist')
    end
else
    switch filter
        case 'LEFT' % <=> RIGHT jacobian
            J11 = sin(th)/th; J12 = (1-cos(th))/th; J13 = (th*v1-v2+v2*cos(th)-v1*sin(th))/th^2;
            J21 = (cos(th)-1)/th; J22 = sin(th)/th; J23 = (v1+th*v2-v1*cos(th)-v2*sin(th))/th^2;
            J = [J11 J12 J13;...
                 J21 J22 J23;...
                  0   0   1];
        case 'RIGHT' % <=> LEFT jacobian
            J11 = sin(th)/th; J12 = (cos(th)-1)/th; J13 = (th*v1+v2-v2*cos(th)-v1*sin(th))/th^2;
            J21 = (1-cos(th))/th; J22 = sin(th)/th; J23 = (-v1+th*v2+v1*cos(th)-v2*sin(th))/th^2;
            J = [J11 J12 J13;...
                 J21 J22 J23;...
                  0   0   1];
        otherwise
            error('filter does not exist')
    end
end

end