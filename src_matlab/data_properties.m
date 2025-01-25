ax_fnt=18;
lab_fnt=18;

sg_change=0
iplot=0
rest_diff=0;
plotseries

maxlead=4;

% y_gr=y_gr(:,ind_ctry);
%ret=find(all((data>-Inf).*(data<Inf),2));

%y_gr=data(ret,:);
y_gr=dlo_unadj;
%cal=cal(ret);
T=size(y_gr,1);
N=size(y_gr,2);
prop_mat=zeros(5,N);
conc_mat=zeros(maxlead+1,N-1,1);

corr_c={'YER'};
c_ind=find(sum(strcmp(ctry(ones(size(corr_c,1),1),:),corr_c(:,ones(1,size(ctry,2)))),1));
corr_l={'KTAUF'};
l_ind=find(sum(strcmp(ctry(ones(size(corr_l,1),1),:),corr_l(:,ones(1,size(ctry,2)))),1));


%y_grma=y_gr-kron(ones(T,1),mean(y_gr,1));
y_grma=dlo;
It=y_grma>=0;
for j=1:N;
    y_graa(j)=mean(y_grma(y_grma(:,j)>=0,j));
    y_grba(j)=mean(y_grma(y_grma(:,j)<0,j));
end
% y=data(:,ind_ctry);

%concordance matrix : (0,..maxlead) x N x (cyclical, countercyclical)
conc_ser={'YER'};
conc_ind=find(sum(strcmp(ctry(ones(size(conc_ser,1),1),:),conc_ser(:,ones(1,size(ctry,2)))),1));
m_It=It(:,[1:conc_ind-1 conc_ind+1:end])-kron(ones(T,1),mean(It(:,[1:conc_ind-1 conc_ind+1:end]),1));
m_It_cycl=kron(It(:,conc_ind)-mean(It(:,conc_ind)),ones(1,size(m_It,2)));
%m_It_anti=kron((1-It(:,conc_ind))-mean((1-It(:,conc_ind))),ones(1,size(m_It,2)));
% m_It_cycl=cat(3,m_It_cycl,m_It_anti);

for l=0:maxlead
    conc_mat(l+1,:,:)=2/(T-l).*sum([m_It(1:end-l,:,:).*m_It_cycl(l+1:end,:,:)],1);
end
%variance of mean adjusted contemporaneous concordance index
[ac, cov, hm, sd, med, ki, ineff, w] = autocovneu_conc(It,10);
cv = conc_var(cov,size(cov,1)-1,conc_ind,w);

if sg_change
    [sg] = neg_corr_sign(y_gr,ctry,corr_c); %change sign for series significantly negatively correlated with GDP
    y_gr=y_gr.*sg(ones(T,1),:);
%     [sg] = neg_corr_sign(y_gr,ctry,corr_l); %change sign for series significantly negatively correlated with leading
%     y_gr=y_gr.*sg(ones(T,1),:);
end
%build matrix of data properties
[R P]=corrcoef(y_gr);
[Rlag Plag]= corrcoef([y_gr(1:end-1,:) y_gr(2:end,:)]); 
% [R8092 P8092]=corrcoef(y_gr(cal(2:end)<1993,:));
% [R9002 P9002]=corrcoef(y_gr(cal(2:end)>=1990,:));
prop_mat=[mean(y_gr,1);y_graa;y_grba;std(y_gr,[],1);R(c_ind,:);P(c_ind,:);R(l_ind,:);P(l_ind,:)];
% mean_mat=[mean(y_gr,1); ... 
%         mean(y_gr(((cal(2:end)>=1980).*(cal(2:end)<1988))==1,:),1);mean(y_gr(((cal(2:end)>=1985).*(cal(2:end)<1993))==1,:),1); ...
%         mean(y_gr(((cal(2:end)>=1990).*(cal(2:end)<1998))==1,:),1);mean(y_gr(((cal(2:end)>=1995).*(cal(2:end)<2003))==1,:),1)];
corr_mat=[R(c_ind,:)' R(l_ind,:)'];
save corr_0703.asc corr_mat -ascii

gh=figure(gcf+1)
scatter(y_grba(or(P(c_ind,:)<0.05,P(l_ind,:)<0.05)),y_graa(or(P(c_ind,:)<0.05,P(l_ind,:)<0.05)),'o');
set(gca,'FontSize',ax_fnt)
xlabel('below average','FontSize',lab_fnt)
ylabel('above average','FontSize',lab_fnt)
print(gh,'-deps',['dat_bamean']);

gh=figure(gcf+1)
scatter(mean(y_gr(:,or(P(c_ind,:)<0.05,P(l_ind,:)<0.05)),1),std(y_gr(:,or(P(c_ind,:)<0.05,P(l_ind,:)<0.05)),[],1),'o');
set(gca,'FontSize',ax_fnt)
xlabel('mean','FontSize',lab_fnt)
ylabel('standard deviation','FontSize',lab_fnt)
print(gh,'-deps',['dat_meanstd']);

gh=figure(gcf+1)
scatter(R(c_ind,or(P(c_ind,:)<0.05,P(l_ind,:)<0.05)),R(l_ind,or(P(c_ind,:)<0.05,P(l_ind,:)<0.05)),'o');
set(gca,'FontSize',16)
xlabel('with GDP','FontSize',ax_fnt)
ylabel('with KTAUF','FontSize',lab_fnt)
print(gh,'-deps',['dat_corr']);


fid=fopen('data_prop.tex','w');
t=['\\begin{table}\\caption{Data properties}\\label{tab:prop}\\begin{center}\\begin{tabular}{l|' ... 
        '*{' int2str(8) '}{c}|} \\hline \\hline \n'];
fprintf(fid,t);

t=['series& mean& mean & mean & stand.dev.& GDP & P-value & KTAUF & P-value \\\\  \n'];
fprintf(fid,t);
t=['& & above av.& below av.& & & & & \\\\ \\hline  \n'];
fprintf(fid,t);

% t2=[];
% t3=[];
for i=1:N;
    t1=[];
    t1=[char(ctry(i)) ];
       for j=1:size(prop_mat,1)-4;
           t1=[ t1 '& ' num2str(prop_mat(j,i),'%5.2f') ];
       end
       if P(1,i)<0.05
           t1=[ t1 '& {\\bf ' num2str(prop_mat(j+1,i),'%5.2f') '}' ];
       else
           t1=[ t1 '& ' num2str(prop_mat(j+1,i),'%5.2f') ];
       end
       t1=[ t1 '& ' num2str(prop_mat(j+2,i),'%5.2f') ];
       if P(l_ind,i)<0.05
           t1=[ t1 '& {\\bf ' num2str(prop_mat(j+3,i),'%5.2f') '}' ];
       else
           t1=[ t1 '& ' num2str(prop_mat(j+3,i),'%5.2f') ];
       end
       t1=[ t1 '& ' num2str(prop_mat(j+4,i),'%5.2f') ];
       t1 = [t1 '\\\\ \n'];
       fprintf(fid,t1);
end
t='\\hline \\hline \\end{tabular}\\end{center}\\end{table}\n \n';
fprintf(fid,t);

t=['\\begin{table}\\caption{Concordance index, mean adjusted (first column, t-statistic of contemporaneous concordance}\\label{tab:conc}\\begin{center}\\begin{tabular}{l|' ... 
        '*{' int2str(maxlead+2) '}{c}|} \\hline \\hline \n'];
fprintf(fid,t);

t=['series&  t-stat. & \\multicolumn{' int2str(maxlead+1) '}{l}{lead} \\\\  \n'];
fprintf(fid,t);
t=['& & 0'];
for j=1:maxlead
    t=[t '& ' int2str(j) ];
end
fprintf(fid,[t ' \\\\ \\hline  \n']);

% t2=[];
% t3=[];
ind_exconc=[1:conc_ind-1 conc_ind+1:N];
for i=1:N-1
    t1=[];
%     t2=[];
    t1=[char(ctry(ind_exconc(i))) '&' num2str(conc_mat(1,i,1)./sqrt(cv(i)),'%5.2f') ];
%     t2=['& a ' ];
    for j=1:maxlead+1;
        t1=[ t1 '& ' num2str(conc_mat(j,i,1),'%5.2f') ];
%         t2=[ t2 '& ' num2str(conc_mat(j,i,2),'%5.2f') ];
    end
    %         t1 = [t1 '\\\\ \n'];
    fprintf(fid,[t1 '\\\\ \n']);
%     fprintf(fid,[t2 '\\\\ \n']);
end
t='\\hline \\hline \\end{tabular}\\end{center}\\end{table}\n \n';
fprintf(fid,t);


fprintf(fid,'P-Value GDP \n');
% t=num2str([P(3,:);P(end,:)],'%8.4f %8.4f\n');
t = [P(1,:)];
fprintf(fid,'%6.4f\n',t);

figure(gcf+1)
sig=sign(conc_mat(find(abs(conc_mat)==kron(ones(maxlead+1,1),max(abs(conc_mat),[],1)))));
scatter(conc_mat(1,:)',sig.*max(abs(conc_mat),[],1)');
xlabel('concordance index (contemporaneous)')
ylabel('maximum concordance index')


% fprintf(fid,'80-02   80-87  85-92  90-97  95-02\n');
% fprintf(fid,'%6.2f  %6.2f %6.2f %6.2f %6.2f\n',mean_mat);

% t=num2str(P(end,:),'%8.4f\n');
% fprintf(fid,'P-Value Corr with US \n');


% figure(gcf+1)
% scatter(R8092(3,[1:2 4:end-1]),R8092(end,[1:2 4:end-1]),'o','filled');
% % scatter(R8092(3,:),R8092(end,:));
% title('correlations of GDP growth 1980-1992')
% xlabel('correlation with German GDP')
% ylabel('correlation with US GDP')
% hold on;
% for j=[1 2 4:N-1];
% % for j=1:N;
%     text(R8092(3,j)-0.01,R8092(end,j),ctries_sel(j,:),'VerticalAlignment','middle','HorizontalAlignment','right','FontSize',12)
% end

% figure(gcf+1)
% scatter(R9002(3,[1:2 4:end-1]),R9002(end,[1:2 4:end-1]),'o','filled');
% % scatter(R8092(3,:),R8092(end,:));
% title('correlations of GDP growth 1990-2002')
% xlabel('correlation with German GDP')
% ylabel('correlation with US GDP')
% hold on;
% for j=[1 2 4:N-1];
% % for j=1:N;
%     text(R9002(3,j)-0.01,R9002(end,j),ctries_sel(j,:),'VerticalAlignment','middle','HorizontalAlignment','right','FontSize',12)
% end

% 

% [corr_sor si]=sort(diag(Rlag(N+1:end,1:N)));
% lab_sor=ctries_sel(si,:);
% figure(gcf+1)
% bar(1:N,corr_sor);
% % scatter(R8092(3,:),R8092(end,:));
% title('correlation of y_{it} with y_{i,t-1}')
% ylabel('correlation')
% hold on;
% for j=[1:N];
% % for j=1:N;
%     text(j,(corr_sor(1)<0)*(corr_sor(1)-0.07),lab_sor(j,:),...
%                'VerticalAlignment','middle','HorizontalAlignment','right','Rotation',90,'FontSize',12)
% end
% 
% t=[ctries_sel num2str(diag(Rlag(N+1:end,1:N)),'%8.4f')  num2str(diag(Plag(N+1:end,1:N)),'%8.4f')]
st=fclose(fid);
% 
% figure(gcf+1)
% for j=[1:N];
% scatter(y_gr(1:end-1,j),y_gr(2:end,j));
% % scatter(R8092(3,:),R8092(end,:));
% hold on;
% end
% title('y_{t-1} against y_t')
% ylabel('y_{t}')
% xlabel('y_{t-1}')
