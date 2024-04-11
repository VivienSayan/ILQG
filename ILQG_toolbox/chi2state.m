function [Rot,theta,x] = chi2state(chi)
Rot = chi(1:2,1:2);
theta = atan2(Rot(2,1),Rot(1,1));
x = chi(1:2,3);
end