function [qq,qmc] = evalq5ddmix2(qqall,qmall,betasim,sepssim,etasim,betamcmc,sepsmcmc,etamcmc,anmc,ancholmc,postseps,etapost)

% Random permutation sampler 
% input qqall, qmall aus der auswertung von q (M times MMC), vergleichbar eval5dd

MMC=size(anmc,2);
L=size(etasim,2);
nst=size(etasim,1);

for l=1:MMC
   
   
      infcholp=inv(squeeze(ancholmc(l,:,:))); %sigma^-1=infchol'*infchol
      
      qmall(:,l) = qmall(:,l) + pmultnormchol(betamcmc,squeeze(anmc(:,l*ones(1,L))),infcholp)';
      qmall(:,l) = qmall(:,l) + squeeze(dirichlog(etamcmc,etapost(l*ones(1,L),:)'))';
      qmall(:,l) = qmall(:,l) + squeeze(pinvgamlog(sepsmcmc,postseps(l*ones(1,L),:)'))';
      
      qqall(:,l) = qqall(:,l) + pmultnormchol(betasim,squeeze(anmc(:,l*ones(1,L))),infcholp)';
      qqall(:,l) = qqall(:,l) + squeeze(dirichlog(etasim,etapost(l*ones(1,L),:)'))';
      qqall(:,l) = qqall(:,l) + squeeze(pinvgamlog(sepssim,postseps(l*ones(1,L),:)'))';
                      
     end

qmax= max(qmall,[],2);
qmc=qmax+log(sum(exp(max(qmall-qmax(:,ones(1,MMC)),-1e15*ones(L,MMC))),2));

qmax= max(qqall,[],2);
qq= qmax+log(sum(exp(max(qqall-qmax(:,ones(1,MMC)),-1e15*ones(L,MMC))),2));

qmc=qmc-log(MMC);
qq=qq-log(MMC);
