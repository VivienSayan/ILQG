function [B] = bcr_ilqg(Hk, Fk, Rk, Qk, P0, N, dimyk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective : Compute the Cramér-Rao lower-bound                          %
%                                                                         %
%               B = bcr(Hk, Fk, Rk, Qk, P0, N, dimyk)                            %
%                                                                         %
% INPUT: Hk = matrice jacobienne de mesure évaluée à la position vraie     %
%        Fk = matrice jacobienne de dynamique évaluée à la position vraie  %
%        Rk = matrice de covariance du bruit de mesure                     %
%        Qk = matrice de covariance du bruit de dynamique                  %
%        P0 = matrice de covariance de l'état initial                      %
%        N = nombres d'itérations  (T/dt)
%        dimyk = dimension du vecteur observé à l'instant k
%                                                                          %
% OUTPUT: B = Cramér-Rao bound selon la formule de Tichavsky 
%
% P. Tichavsky, C. H. Muravchik and A. Nehorai, "Posterior 
% Cramer-Rao bounds for discrete-time nonlinear filtering,"        
% in IEEE Transactions on Signal Processing, vol. 46, no. 5, pp. 1386-1396
%, May 1998    
%
% Karim Dahia, IGNC                                                       %
% October 2024                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation
B = zeros(length(P0),length(P0),N);
B(:,:,1) = P0;
jP = inv(B(:,:,1));  
    
% Calcul de la BCR
for t = 2:N
    dim_obs = dimyk(t);
    F = Fk(:,:,t);
    H = Hk(1:dim_obs,:,t);
    R = Rk(1:dim_obs,1:dim_obs,t);
    Q = Qk(:,:,t);
    % F = Fk;
    % H = Hk;
    % R = Rk;
    % Q = Qk;
    D11 = F'*inv(Q)*F; % sans l'espérance mathématique
    D12 = -F'*inv(Q);
    D21 = D12';
    D22 = inv(Q) + H'*inv(R)*H;
    jP = D22 - D21*inv(jP + D11)*D12;
    B(:,:,t) = inv(jP);
end
end
    
  

