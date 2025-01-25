function fd = pinvgamlog(theta,par)

%input: theta ... argument (vector) 
%       par ... parameter (zeilenvector)


alpha = par(1,:);
beta = par(2,:);
fd=alpha.*log(beta)-gammaln(alpha)-beta./theta-(alpha+1).*log(theta);
