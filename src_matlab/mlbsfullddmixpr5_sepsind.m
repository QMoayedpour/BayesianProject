function [mllogbs, sdrel,loglikq,loglikmc,priorq,priormc,qq,qmc,mu,Qsim,Qsiminv,detQsiminv,qmall,qqall] = mlbsfullddmixpr5_sepsind(alphamc,sepsmc,etamc,Qmc,Qinvmc,Qinvdet,prQnu,prQS,postQnu,postQS,e0,etapost,prs,postseps,prlambda,anmc,ancholmc,pralm,pralinf,y,ZG,dd)
% bridge sampling estimators: random permuations  from MCMC sample of size MMC (number of permuations: see evalq4)


%NEU: prQnu, prQS, Qmc, postQnu,postQS, ÜBERGEBEN

%r=size(W,2);
%prQS=matrix r r;
%prQnu scalar

M=size(sepsmc,1); %  laenge MCMC sample
L=M;      % laenge q sample
MMC=size(postseps,1);  % laenge der mixture
% MMC=10;  % zum testen
% joint prior beta

% Call = C0;
% call = c0;
% bd = size(c0,1);
bd = 0;
r = size(ZG,2);
nst = size(etamc,2);


% je nach prior aufblasen oder kommentieren:

d=r;
K=nst;

% index1 = [1:d];
% index = reshape(index1(ones(1,K),:)',1,K*d);
% states = [1:K]';
% index2 = reshape(states(:,ones(1,d))',1,K*d);

index1 = [1:dd];
index = reshape(index1(ones(1,K),:)',1,K*dd);
states = [1:K]';
index2 = reshape(states(:,ones(1,dd))',1,K*dd);

pralm = pralm([index dd+1:d],1);                            
prcov = zeros(d+(K-1)*dd,d+(K-1)*dd);
for k=1:K,
    prcov(1+(k-1)*dd:k*dd,1+(k-1)*dd:k*dd) = pralinf(1:dd,1:dd);
end 
if d > dd
    prcov(K*dd+1:d+(K-1)*dd,K*dd+1:d+(K-1)*dd)=pralinf(dd+1:d,dd+1:d);   % fuer X1==X2
end   
pralinf=prcov;
%%%%%%%%%%%%%%%%%%%

e0=e0';

r = size(ZG,2);
nst = size(etamc,2);
% kfac = gammaln(nst+1);
ancholsize=size(ancholmc,2);
betasim=zeros(ancholsize,L);
sQ=size(Qmc);
Qsim=zeros(size(Qmc));
Qsiminv=zeros(size(Qmc));
detQsiminv=zeros(size(squeeze(Qmc(:,1,:))));




% for i = 1:nst
%    call = [call; prm(:,i)];
%    Fillzeros = zeros(bd+(i-1)*r,r);
%    Call = [Call Fillzeros; Fillzeros' prv(:,:,i)];
% end

% 1a. construction of q and sampling from q

%    (Wahl von q: mixture approximation)


anmc=anmc';

% 1b. simulate from q 

'1.a: simulate q'

mu = simuni(MMC,L);  
sepssim = 1./gamrnd(squeeze(postseps(mu,1)),1./squeeze(postseps(mu,2)))';

gam = gamrnd(etapost(mu,:)',ones(nst,L));   
sumgam =  sum(gam,1);
etasim =  gam./ sumgam(ones(1,nst),:);


%scale factor for conditional  density 
comp=[1:M]';

if any(any(squeeze(Qmc(1,:,:))))
    detprQS=log(det(prQS));
else
    detprQS=0;
end   
prQScol=qincol(prQS);
detS=zeros(MMC,K);


for m=1:MMC
    %   '1.a'
    m
    index=mu==m;
    lm = sum(mu==m);
    cc=comp(index);
    eps = randn(ancholsize,lm);
    betasim(:,mu==m) = squeeze(anmc(:,ones(1,lm)*m)) + squeeze(ancholmc(m,:,:))*eps;
    
    if any(any(squeeze(Qmc(1,:,:))))
        
        for k=1:nst 
            S=qinmatr(postQS(m,:,k)');
            [us t] = schur(S);
            
            Tdiag=max(diag(real(t)),0);
            detS(m,k) = sum(log(Tdiag));
            
            th=diag(Tdiag.^.5);
            u=real(us)*th;
            
            thinv=diag(Tdiag.^(-.5));
            uinv=real(us)*thinv;
            
            for l=1:lm
                % Qsim(cc(l),:,k) = qincol(u*raninvwi_eye(postQnu(m,k),size(th,1))*u');
                % simuliere inverse
                E=ranwi_eye(postQnu(m,k),size(thinv,1));
                Qsim(cc(l),:,k) = qincol(u*inv(E)*u');
                Qsiminv(cc(l),:,k) = qincol(uinv*E*uinv');
                detQsiminv(cc(l),k)=log(det(E))-detS(m,k);
            end   
            
        end
    end
end   

etamcmc=etamc';
sepsmcmc=sepsmc';
alphamcmc=alphamc';

% 1c. evaluate prior and marglik for q sample 


loglikq=zeros(L,1);
loglikmc=zeros(L,1);
priorq=zeros(L,1);
priormc=zeros(L,1);

priorq = squeeze(dirichlog(etasim,e0(:,ones(1,L))))';
priorq = priorq+pmultnormlog(betasim,pralm,pralinf)';

priormc = squeeze(dirichlog(etamcmc,e0(:,ones(1,L))))';
priormc = priormc + pmultnormlog(alphamcmc,pralm,pralinf)';

if all(prs~=0)
 priorq = priorq+squeeze(pinvgamlog(sepssim,prs(ones(1,L),:)'))';  
 priormc = priormc+squeeze(pinvgamlog(sepsmc',prs(ones(1,L),:)'))'; 
elseif  any(prs~=0)
    'use of improper prior on seps not possible when computing model likelihoods' 
end


qqall=zeros(L,MMC);
qmall=zeros(L,MMC);


Qm0=qinmatr(Qsim(1,:,1)');
Qm=zeros(size(Qm0,1),size(Qm0,2),nst);
Qq=Qm;
Qqlik=Qq;
Qmlik=Qm;

'1.b: evaluate likelihood'


for m=1:L
    %  '1.b: ' 
    if any([1:1000:L]==m)
        m
    end
    if any(any(squeeze(Qmc(1,:,:))))
        
        
        for k=1:nst
            
            Qqlik(:,:,k)=qinmatr(Qsim(m,:,k)');
            Qq1=squeeze(Qsiminv(m,:,k))';
            detQq1 = detQsiminv(m,k);
            
            priorq(m) = priorq(m) + pinvwilog_neu(Qq1,detQq1,prQnu,prQScol,detprQS);
            
            
            Qmlik(:,:,k)=qinmatr(Qmc(m,:,k)');
            Qm1=squeeze(Qinvmc(m,:,k))';
            detQm1 = Qinvdet(m,k);
            
            
            priormc(m) = priormc(m) + pinvwilog_neu(Qm1,detQm1,prQnu,prQScol,detprQS);
            
            
            qmall(m,:) = qmall(m,:)+ pinvwilog_neu(Qm1(:,ones(1,MMC)),detQm1,postQnu(:,k)',squeeze(postQS(:,:,k))',detS(:,k)');
            qqall(m,:) = qqall(m,:)+ pinvwilog_neu(Qq1(:,ones(1,MMC)),detQq1,postQnu(:,k)',squeeze(postQS(:,:,k))',detS(:,k)');
            
        end
    end
    
    loglikq(m) = simstateexstudent(y,etasim(:,m)',ZG,Qqlik,betasim(:,m)',sepssim(1,m),prlambda,dd);
    loglikmc(m) = simstateexstudent(y,etamcmc(:,m)',ZG,Qmlik,alphamcmc(:,m)',sepsmcmc(m),prlambda,dd);
    
end   

[qq,qmc] = evalq5ddmix2(qqall,qmall,betasim,sepssim,etasim,alphamcmc,sepsmcmc,etamcmc,anmc,ancholmc,postseps,etapost);

ratio = loglikq+priorq-qq;
ratiomax=max(ratio);
ratioexp=exp(ratio-ratiomax);

mllogmi = ratiomax+log(mean(ratioexp)); 
mlsdmi = (L)^(-.5)*exp(log(std(ratioexp))-log(mean(ratioexp)));

sdrel(1)= mlsdmi;

% 1d. berechnen der IS - modellikelihood

% 2a. Evaluate q at the MCMC sample

%etamcmc=etamc(it0+1:it,:)';
%sepsmcmc=sepsmc(it0+1:it)';
%alphamcmc=alphamc(it0+1:it,:)';

% priormc=zeros(L,1);


% 2b. Berechnen der RIS- modelliekihood


r=qmc-loglikmc-priormc;
[ls is]=sort(r);
logc=-r(is(end));

rqp=exp(max(r-max(r),-1e15))';
mllogri = -max(r)-log(mean(rqp));

[omega ac eib abserrn eff] = autocov_alt(rqp,1000);
mlsdri = (abserrn(1,1))^.5/mean(rqp);
sdrel(2)=mlsdri;

% 3. iterative bridge  sampling, starting from ri and mi

maxit=50;
mllogbs=zeros(maxit,2);
mllogbs(1,1)=mllogmi;  
mllogbs(1,2)=mllogri;  

for i=2:maxit
    %  i  
    logpostmc=loglikmc+priormc-mllogbs(i-1,2);
    logpostq=loglikq+priorq-mllogbs(i-1,2);
    maxqmc=max([qq logpostq],[],2);
    rq=logpostq-maxqmc-log(exp(log(M)+logpostq-maxqmc)+exp(log(L)+qq-maxqmc));
    maxqmc=max([qmc logpostmc],[],2);
    rmc=qmc-maxqmc-log(exp(log(M)+logpostmc-maxqmc)+exp(log(L)+qmc-maxqmc));
    %mllogbs(i,2)=mllogbs(i-1,2)+log(mean(rq)/mean(rmc));
    mllogbs(i,2)=mllogbs(i-1,2)+max(rq)+log(mean(exp(max(rq-max(rq),-1e15))))-max(rmc)-log(mean(exp(max(rmc-max(rmc),-1e15))));
    
    
    logpostmc=loglikmc+priormc-mllogbs(i-1,1);
    logpostq=loglikq+priorq-mllogbs(i-1,1);
    maxqmc=max([qq logpostq],[],2);
    rq=logpostq-maxqmc-log(exp(log(M)+logpostq-maxqmc)+exp(log(L)+qq-maxqmc));
    maxqmc=max([qmc logpostmc],[],2);
    rmc=qmc-maxqmc-log(exp(log(M)+logpostmc-maxqmc)+exp(log(L)+qmc-maxqmc));
    mllogbs(i,1)=mllogbs(i-1,1)+max(rq)+log(mean(exp(max(rq-max(rq),-1e15))))-max(rmc)-log(mean(exp(max(rmc-max(rmc),-1e15))));
    
    
    
end   

%standard deviations bridge sampling:

rmc=exp(max(rmc-max(rmc),-1e15));
rq=exp(max(rq-max(rq),-1e15));

[omega ac eib abserrn eff] = autocov_alt(rmc',1000);
mlsdbs =(abserrn(1,1)/(mean(rmc)^2)+((std(rq)/mean(rq))^2)/L)^.5;
sdrel(3)=mlsdbs;
