ifig=ifig+1;figure(ifig);
if var_logit==1
load canova_K2_end1_prlambda8_Slogit1.mat;
end
z1=min(Zlogit(:,2));
z1=min(z1,-1.2);
z2=max(Zlogit(:,2));
zz=linspace(z1,z2);
gamz=gammc(:,2,2);
pr_logit=zeros(size(zz));
for i=1:size(zz,2)
ee=gammc(:,1,2)+gammc(:,2,2)*zz(1,i);
pr_logit(1,i)=mean(1./(1+exp(ee)));
end
plot(zz,pr_logit);
hold on;
scatter(Zlogit(:,2),zeros(size(Zlogit(:,2))),3,'filled');
print  -depsc2 can_class_prior_logit.eps;

pr_logit=zeros(size(Zlogit(:,2)));
for i=1:size(Zlogit(:,2),1)
    ee=gammc(:,1,2)+gammc(:,2,2)*Zlogit(i,2);
pr_logit(i)=mean(1./(1+exp(ee)));
end

ifig=ifig+1;figure(ifig);
load canova_K2_end1_prlambda8.mat;
pr=mean(etamc(:,1));
plot([z1 z2],pr*ones(1,2),'-');
print  -depsc2 can_class_prior.eps;


ifig=ifig+1;figure(ifig);
if var_logit==1
load canova_K2_end1_prlambda8_Slogit1.mat;
end
run_model_endswit;
pp=zeros(N,1);
pp_logit=zeros(N,1);
M=size(lambdamc,1);
MM=size(alphamc,1);dM=MM-M;
   for i=1:N,
    for m=dM+1:MM
       for k=1:K, Q0(:,:,k) = qinmatr(Qmc(m,:,k));end
        pmi=prob_mix_logit(y(:,i),squeeze(gammc(m,:,:)),Z(:,:,i),Zlogit(i,:),Q0,alphamc(m,:),sepsmc(m,:)/lambdamc(m-dM,i),dd,0);  
        pp(i) = pp(i)+pmi(1);  
  end
end
pp_logit=pp/M;
 scatter(Zlogit(:,2),pp_logit,3,'filled');
% plot(Zlogit(:,2),pp_logit);
 pp=zeros(N,1);
 load canova_K2_end1_prlambda8.mat;
M=size(lambdamc,1);
MM=size(alphamc,1);dM=MM-M;
   for i=1:N,
    for m=dM+1:MM
       for k=1:K, Q0(:,:,k) = qinmatr(Qmc(m,:,k));end
        pmi=prob_mix(y(:,i),squeeze(etamc(m,:)),Z(:,:,i),Q0,alphamc(m,:),sepsmc(m,:)/lambdamc(m-dM,i),dd,0);  
        pp(i) = pp(i)+pmi(1);  
  end
end
pp=pp/M;
S_logit=1;run_model_endswit;
 scatter(Zlogit(:,2),pp,3,'filled','d');
 %plot(Zlogit(:,2),pp,'.--');

print  -depsc2 can_class_pos.eps;

tab=[pr_logit pp_logit  pp];
maketab(tab);
save canova_K2_discrimate pr_logit pp_logit pr pp;
