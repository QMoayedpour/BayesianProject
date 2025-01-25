function [omega, ac, hm, varhhat, eff] = autocov(h,r)
%
% h ... anzahl der zeilen = anzahl der variablen
% r ... maximaler lag
% output:
% omega(k,k,1:r+1)
%   k ... variablen
%   omega(:,:,1) ... covariance matrix
%   omega(:,:,s) ... autocov matrix for lag s-1
% ac(1:r,k) 
%   k .. anzahl der variablen
%   ac(s,k) autocorrelation at lag s
% hm(1:k)
%   hhat(j) ... MW der j.ten variable
% varhhat(1:k,1:k)
% eff(1:k) .. efficiency compared to i.i.d mean
%   
% um die ineff richtig vergleichen zu koennen, muessen die ac in das Band 
% hineinkommen!!! KONTROLLIEREN !!!!

n = size(h,2);
k = size(h,1);
hm = mean(h,2);
omega = zeros(k,k,r+1);
omegat = omega;
svec = zeros(k,k,r+1);
ac = zeros(r,k);
eff = zeros(k,1);

for j=1:k
hhat(j,1:n) = hm(j);
end

for s = 1:r+1
omegas(1:k,1:k) = (h(:,s:n)-hhat(:,1:n-s+1))*(h(:,1:n-s+1)-hhat(:,1:n-s+1))'./n;
omega(:,:,s) = omegas(:,:);
omegat(:,:,s) = omegas(:,:)';
svec(1:k,1:k,s) = (1-(s-1)/(r+1))/n;
end
svec(:,:,1) = svec(:,:,1)*0.5;

for j=1:k
ac(:,j) = squeeze(omega(j,j,2:end))./omega(j,j,1);
end

varhhat = squeeze(sum(svec.*(omega+omegat),3));
for j=1:k
eff(j) = varhhat(j,j)/omega(j,j,1)*n;
end

