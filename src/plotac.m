function ineff= plotac(x,r,cl)

% x ... observations
% r ... max lag of ac
% cl ... color of plot

[ac, hm, sd, med, ki, ineff] = autocovneu(x,r);
plot([0:r]',[1;ac],cl);
axis([0 r -0.05 1.05])
hold on
sd=2/sqrt(size(x,1));
plot([0 r],-sd*[1 1],'k');
plot([0 r],sd*[1 1],'k');
