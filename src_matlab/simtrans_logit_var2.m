function  [gam,mr,sr] = simtrans_logit_var2(Zlogit,gam,prgamm,prgaminf,D,mr,sr);
  
K = size(gam,2);
dg=size(gam,1);

% 1. Sample hidden propensities based on gam old

lamlog=Zlogit*gam;
lam=exp(lamlog); % predictor

 U=rand(size(lam)).*(1-D')+D'.*ones(size(lam));
 ymin= -log(rand(size(lam,1),1))./sum(lam,2);
ystar= ymin(:,ones(1,K))- log(U)./lam;
ystarlog=log(ystar);
ym=(ystarlog-mr)./sr;

% 2. Sample new parameter of the logit model

index=[1:K];
 k0=index(all(gam==0,1)==1);k0=k0(1);   % baseline category k0 

baseline=0; %baseline=0: sample also baseline category from unidentified model
baseline=1; %baseline=1: fix baseline category k0 (no permutation sampling : 1);

if baseline
    indexk=[1:k0-1 k0+1:K];
else
    indexk=index;
end


for k=indexk
    Zk=-Zlogit./sr(:,k*ones(1,dg));
    postgaminf=prgaminf+Zk'*Zk;
    postgamcov=inv(postgaminf);
    postgamm=postgamcov*(Zk'*ym(:,k)+prgaminf*prgamm);
    gam(:,k)=postgamm+chol( postgamcov)'*randn(dg,1);
end    


% 3. Sample component indicator m and S

gam=gam-gam(:,k0*ones(1,size(gam,2)));

 [mr,sr,r]=draw_comps_mix(ystarlog+Zlogit*gam);

 