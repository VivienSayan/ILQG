function res = dl_deps(eps)
eps_th = eps(1); eps_x = eps(2); eps_y = eps(3);

dl_deps_th = (eps_x*cos(eps_th)-eps_y*sin(eps_th))/eps_th - (eps_y*cos(eps_th)-eps_y+eps_x*sin(eps_th))/eps_th^2;
dl_deps_x = sin(eps_th)/eps_th;
dl_deps_y = (cos(eps_th)-1)/eps_th;

res = [dl_deps_th, dl_deps_x, dl_deps_y];

end