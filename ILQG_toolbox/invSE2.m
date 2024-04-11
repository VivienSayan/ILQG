function chi_inv = invSE2(chi)

R = chi(1:2,1:2);
t = chi(1:2,3);
v = -R'*t;

chi_inv = [[R',v];[zeros(1,2),1]];

end