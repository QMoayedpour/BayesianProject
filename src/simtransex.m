function [eta,etapost] = simtransex(D,c0)

K = size(D,1);

ct = sum(D,2)';
etapost = ct + c0;
gam = gamrnd(etapost,1);
eta = gam ./ (sum(gam,2));
