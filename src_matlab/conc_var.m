function cv = conc_var(cov,tau,conc_ind,w);

% function to estimate the variance of an adjusted concordance index
% see Artis et al. (2004), Dating business cycles, Oxford Bulletin of
% Economics & Statistics, 66, 537-565.
% cov .. matrix of covariance functions (tau+1) x N(umber of series)
% tau .. truncation parameter
% conc_ind .. index of the series relative to which the concordance index
%             is computed
% sigma^2 = gamma(0)_i  gamma(0)_j + 2 \sum (1-tau/T) gamma(tau)_i  gamma(tau)_j

N=size(cov,2)-conc_ind;
%w = kron([1 2*(1-[1:tau]./(tau+1))]',ones(1,N));
w=2*w(:,ones(1,N));
    
cv = sum((kron(cov(:,conc_ind),ones(1,N)).*w) .* cov(:,[conc_ind+1:end]),1).*4;