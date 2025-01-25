
load dating 
pp=real(pp);
ecri=real(ecri);
pds=real(pds);
bds=real(bds);
uni=real(uni);

fore_ser=mean(bfore_mat,3);
calp=cal(end-size(dlo,1)+1:end);
calf=calp(1):1/freq:cal_fore(end);
cald=cal_dat(1):1/freq:max(cal_fore(end),cal_dat(end));
ff='%5.2f ';
for h=1:9;
    ff=[ff '%5.2f '];
end
ff=[ff '\n'];

for g=1:size(bfore_mat,1)
    figure%(gcf+1)
    for h=1:H+1
    plot(cal_fore(h:h+Hf-1),fore_ser(g,:,:,h),'r')
    hold on
    end
    plot(calp,dlo(:,vars(g)))
    axis([calf(45) calf(end) min(dlo(:,vars(g)))-0.05  max(dlo(:,vars(g)))+0.05]);
%    legend('2004.75','2005.0','2005.25','2005.5','2005.75','actual')
    title(ser_fc(g))
end


probMS=1-mean(LIMSmc,3);
prob_fore=1-mean(bforeIMS_mat,3);
calp=cal(end-size(probMS,2)+1:end);
calf=calp(1):1/freq:cal_fore(end);
%for g=1:ceil(size(bforeIMS_mat,4)/4)
figure%(gcf+1)
for h=1:H+1
    subplot(1+(K>1),1,1)
    bar(calp,probMS(1,:),'c')
    hold on
    plot(calf,ones(length(calf),1)*0.5);
    hold on
%     plot(cal_fore(h:h+Hf-1),squeeze(prob_fore(1,:,:,h)),'r')
    plot(cal_fore(h+1-1),squeeze(prob_fore(1,1,:,h)),'r+')
    hold on
    axis([calf(1) calf(end) 0 1])
    title(['forecast horizon: ' num2str(cal_fore(h),'%5.2f') '--' num2str(cal_fore(h+Hf-1),'%5.2f')])
    if K>1
    subplot(2,1,2)
    bar(calp,probMS(2,:),'c')
    hold on
    plot(calf,ones(length(calf),1)*0.5);
    hold on
%    plot(cal_fore(h:h+Hf-1),squeeze(prob_fore(2,:,:,h)),'r')
    plot(cal_fore(h+1-1),squeeze(prob_fore(2,1,:,h)),'r+')
    hold on
    axis([calf(1) calf(end) 0 1])
    end
end
figure%(gcf-(g+1))
probMS(:,end+1:end+Hf)=0/0;
c_adj=length(cald)-length(cal_dat);
ecri(end+1:end+c_adj)=0/0;
pds(end+1:end+c_adj)=0/0;
pp(end+1:end+c_adj)=0/0;
bds(end+1:end+c_adj)=0/0;
uni(end+1:end+c_adj)=0/0;

% t=['forecast date; YER and KTAUF growth; prob_fore (2); prob_insample(2); in-sample comparison(2); versus ECRI(1) ; versus PDS(1) \n'];
t=['forecast date; prob_fore (1); prob_insample(1); in-sample comparison(1); versus ECRI(1) ; ECRI; versus PDS(1) ;PDS ; versus PDS(1) ;PDS\n'];

fprintf(fid,t);

for j=1:size(fore_ser,4)
    cal_ind=find(any(calf(ones(Hf,1),:)' == cal_fore(j:j+Hf-1,ones(length(calf),1))',2));
    cal_dind=find(any(cald(ones(Hf,1),:)' == cal_fore(j:j+Hf-1,ones(length(cald),1))',2));
%     t=[cal_fore(j:j+Hf-1)';fore_ser(:,:,:,j);prob_fore(:,:,:,j);probMS(1:2,cal_ind);...
%         (prob_fore(:,:,:,j)>0.5)==(probMS(1:2,cal_ind)>0.5);(prob_fore(1,:,:,j)>0.5)==(ecri(cal_ind)>0.5)';...
%          (prob_fore(1,:,:,j)>0.5)==(pds(cal_ind)>0.5)'];
    t=[cal_fore(j:j+Hf-1)';prob_fore(1,:,:,j);probMS(1,cal_ind);...
        (prob_fore(1,:,:,j)>0.5)==(probMS(1,cal_ind)>0.5);(prob_fore(1,:,:,j)>0.5)==(ecri(cal_dind)>0.5)';
         ecri(cal_dind)' ; ...
         (prob_fore(1,:,:,j)>0.5)==(pds(cal_dind)>0.5)'; pds(cal_dind)';...
         (prob_fore(1,:,:,j)>0.5)==(pds(cal_dind)>0.5)'; pds(cal_dind)'];
fprintf(fid,ff,t);
fprintf(fid,' \n');
end
