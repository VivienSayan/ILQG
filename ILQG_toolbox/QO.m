function [NewPts,Weights] = QO(mu,kmax,Xmean,P,Pts,dim_tab)
% Optimal Quantization algorithm (Pdf = normal distribution) :
% INPUTS :
    % - mu   : learning rate
    % - kmax : maximum number of iteration
    % - Xmean: mean of the distribution
    % - Px   : covariance matrix of the distribution                                                        
    % - Pts  : Initial Sigma-Points Xj, j=1,2*nx+1
    % - dim_tab : index of dimension(s) to update
% OUTPUTS :
    % - NewPts : Updated Points

NewPts = Pts - Xmean;
random_tab = zeros(length(dim_tab),1) + sqrtm(P(dim_tab,dim_tab))*randn(length(dim_tab),kmax);
counter = zeros(1,size(Pts,2));
start_count = 0;
for k = 1:kmax
    Xrand = random_tab(:,k);
    Xrand_tab = repmat(Xrand,1,size(NewPts,2));
    norm2 = vecnorm(Xrand_tab-NewPts(dim_tab,:),2,1);
    [~,winning_index] = min(norm2);
    NewPts(dim_tab,winning_index) = NewPts(dim_tab,winning_index) - abs(mu/k)*(NewPts(dim_tab,winning_index)-Xrand);
    if k >= start_count
        counter(winning_index) = counter(winning_index)+1;
    end
end
NewPts = Xmean + NewPts;
Weights = counter./(kmax-start_count);


end