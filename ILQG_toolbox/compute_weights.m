function [Wm,Wc,ksi,KSI,lambda] = compute_weights(n,alpha,beta,kappa)

lambda = alpha^2*(n+kappa)-n;
Wm = 1/(2*(n+lambda))*ones(1,2*n+1); Wc = Wm;
Wm(1) = lambda/(n+lambda); Wc(1) = Wm(1) + (1-alpha^2+beta);
ksi = sqrt(n+lambda);
KSI = diag(ksi*ones(n,1));

end