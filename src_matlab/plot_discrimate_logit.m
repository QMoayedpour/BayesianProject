load canova_K2_end1_prlambda8_Slogit.mat;
z1=min(Zlogit(:,2));
z1=min(z1,-1.2);
z2=max(Zlogit(:,2));
zz=linspace(z1,z2);
gamz=gammc(:,2,2);
prob=zeros(size(zz));
for i=1:size(zz,2)
ee=gammc(:,1,2)+gammc(:,2,2)*zz(1,i);
prob(1,i)=mean(1./(1+exp(ee)));
end
plot(zz,prob);
hold on;
scatter(Zlogit(:,2),zeros(size(Zlogit(:,2))),3,'filled')

load canova_K2_end1_prlambda8.mat;
ee=mean(etamc(:,1));
plot([z1 z2],ee*ones(1,2),'--')
