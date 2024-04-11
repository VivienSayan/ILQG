function y = h(chi,r) 
% fonction d'observation h() repris depuis le papier suivant:
% "Unscented Kalman Filtering on Lie Groups (M.Brossard, S.Bonnabel, J-P Condomines, IROS 2017)"

y = chi*[0;0;1] + [r;0];
y = y(1:2);

end