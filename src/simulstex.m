function S = simulstprior(N,eta);

% simuliere states aus der prior
% function call

rnd = rand(N,1);
S = zeros(1,N);

for i = (N):-1:1
    p = eta;
    st = sum(cumsum(p) < rnd(i)) + 1;      % sampling
    S(i) = st;
end


