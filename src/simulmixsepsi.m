% simulmix
% load simmix2lcmdat0
load Z
alpha=[4 5 3];
eta=[1];
seps=250;nu=30;
Q=diag([5 5 5]);

N=size(Z,3);
T=size(Z,1);
K=size(eta,2);
dd=size(Q,1);
d=size(Z,2);

alphai=zeros(d,K);

% simuliere die states
S = simulstex(N,eta);

% simuliere seps(i);

lambda=gamrnd(nu/2,ones(N,1)*2/nu)';
sepsi=seps(:,ones(1,N))./lambda;

y = sepsi(ones(1,T),:).^.5.*randn(T,N);
alphai(1:dd,1:K) = reshape(alpha(1,1:dd*K),dd,K);
a=alpha(1,dd*K+1:d+(K-1)*dd)';
alphai(dd+1:d,:) = a(:,ones(1,K));

for i=1:N,
   
   % simuliere die random-effects
   b(:,i) = chol(Q)*randn(dd,1);
   
   % simuliere die Daten
   y(:,i) = y(:,i) + squeeze(Z(:,:,i))*alphai(:,S(i))+squeeze(Z(:,1:dd,i))* b(:,i);
end   