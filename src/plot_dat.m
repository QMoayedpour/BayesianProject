load dating0802
sn=6;

figure(gcf+1)
subplot(sn,1,1)
bar(cal_dat,ecri,0.2)
ylabel('ECRI')
axis([cal_dat(1)-0.25 cal_dat(end)+0.25 0 1])
subplot(sn,1,2)
bar(cal_dat,pds,0.2)
ylabel('PDS')
axis([cal_dat(1)-0.25 cal_dat(end)+0.25 0 1])
subplot(sn,1,3)
bar(cal_dat,itl,0.2)
axis([cal_dat(1)-0.25 cal_dat(end)+0.25 0 1])
ylabel('ITL')
subplot(sn,1,4)
bar(cal_dat,pp,0.2)
ylabel('PP')
axis([cal_dat(1)-0.25 cal_dat(end)+0.25 0 1])
subplot(sn,1,5)
bar(cal_dat,bds,0.2)
ylabel('BDS')
axis([cal_dat(1)-0.25 cal_dat(end)+0.25 0 1])
subplot(sn,1,6)
bar(cal_dat,uni,0.2)
ylabel('UNI')
axis([cal_dat(1)-0.25 cal_dat(end)+0.25 0 1])




