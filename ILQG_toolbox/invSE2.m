function chi_inv = invSE2(chi)

R   = chi(1:2,1:2);
tmp = chi(1:2,3);
v = -R'*tmp;

chi_inv = [[R',v];[zeros(1,2),1]];

end