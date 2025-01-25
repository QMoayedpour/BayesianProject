function fl = dirichlog(eta,par)

% input
% eta ... argument der aposteriori dichte
% par ... parameter der aposteriori dichte

fl=gammaln(sum(par,1))-sum(gammaln(par)-(par-1).*log(eta),1);

