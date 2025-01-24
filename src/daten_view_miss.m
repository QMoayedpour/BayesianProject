%iplot=0;ifig=20;nplot=12;nline=3;xmin=-50;xmax=100;
iplot=0;ifig=20;nplot=4;nline=2;xmin=-5;xmax=5;
nobs=size(dlo_obs,1);
cend = 2004; %with DAT_roh_Q0410.xls
% cend = 2006.25;
%with DAT_roh_Q0412.xls
cal = [sort(cend-[1:nobs]/4)]';

Tall=cal(end-size(dlo,1)+1:end);
dlo_max=max(dlo_obs);dlo_min=min(dlo_obs);
index=[1:nbank];
imiss=zeros(size(index));
% ctries=['AUT';'BEL';'DEU';'ESP';'FIN';'FRA';'GRC';'IRL';'ITA';'LUX';'NLD';'PRT';'DNK';'GBR';'SWE';'CHE';'NOR';'AUS';'CAN';'JPN';'USA'];

N=size(indexb,2);
      for i=1:N
         ib=indexb(i);
         imiss(ib)=any(dlo_out(:,ib)'==1);
     end
imiss_ind=find(imiss==1);
%var=[];
for i=index(imiss(1:imiss_ind(end))==1)
   iplot=iplot+1;
   if iplot==1;ifig=ifig+1;figure(ifig);end
   subplot(nplot/nline,nline,iplot);
   plot(Tall,dlo_obs(:,i)-(dlo_obs(1,i)-dlo_miss(1,i)),'k--',Tall,dlo_miss(:,i),'k');
%   var=[var dlo_obs(:,i) dlo_miss(:,i)];
   axis([Tall(1) Tall(end)  xmin xmax]);
   alpha=0.3;
   %axis([0 nobs+2  dlo_min(i)/(1-alpha) dlo_max(i)*(1+alpha)]);
%    legend(ctries(i,:));
   if iplot==nplot;iplot=0;end
end   
