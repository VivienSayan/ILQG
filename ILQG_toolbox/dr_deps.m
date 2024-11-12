function res = dr_deps(eps)
eps_th = eps(1); eps_x = eps(2); eps_y = eps(3);

dr_deps_th = (eps_y*cos(eps_th)+eps_x*sin(eps_th))/eps_th - (eps_x-eps_x*cos(eps_th)+eps_y*sin(eps_th))/eps_th^2;
dr_deps_x = (1-cos(eps_th))/eps_th;
dr_deps_y = sin(eps_th)/eps_th;

res = [dr_deps_th, dr_deps_x, dr_deps_y];

end