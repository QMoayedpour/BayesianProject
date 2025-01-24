fid=fopen('conc_ind.tex','w');

maxlead=4;
ctry={'ECRI_I';'ECRI_II';'PanelDS';'ITL';'PooledPanel'; 'BivarDS';'Univar'};
I_end=2006.75;

load dating
cal=cal_dat;

IND=1;
if IND==1;
    TP_mat=[1990.0 1993.0; 1994.75 1996.0; 1998.25 1999; 1999.5 2001.75;2003 2003.75; 2004.5 cal(end)];
elseif IND==2;
    TP_mat=[1990.0 1993.25;1994.75 1996.25;1998.25 1999; 1999.25 2002;2002.75 2004;2004.5 cal(end)];
end

E_I=zeros(size(cal));
for j=1:size(TP_mat,1);
    E_I=E_I+and(cal>TP_mat(j,1),cal<=TP_mat(j,2));
end
% load dating.asc %first ECRI dating, second Panel DS
It=[ecri E_I pds itl pp bds uni];
N=size(It,2);
It=It(cal<=I_end,:);
T=length(It);

%concordance matrix : (0,..maxlead) x N x (cyclical, countercyclical)
for conc_ind=1:3;
    m_It=It(:,[conc_ind+1:end])-kron(ones(T,1),mean(It(:,[conc_ind+1:end]),1));
    m_It_cycl=kron(It(:,conc_ind)-mean(It(:,conc_ind)),ones(1,size(m_It,2)));
    %m_It_anti=kron((1-It(:,conc_ind))-mean((1-It(:,conc_ind))),ones(1,size(m_It,2)));
    % m_It_cycl=cat(3,m_It_cycl,m_It_anti);

    conc_mat=zeros(maxlead+1,N-conc_ind,2);

    for l=0:maxlead
        conc_mat(l+1,:,1)=(2/(T-l)).*sum([m_It(1:end-l,:,:).*m_It_cycl(l+1:end,:,:)],1);
        conc_mat(l+1,:,2)=(1/(T-l)).*sum([It(1:end-l,conc_ind*ones(1,N-conc_ind),:).*It(l+1:end,conc_ind+1:end,:)]+...
            [(1-It(1:end-l,conc_ind*ones(1,N-conc_ind),:)).*(1-It(l+1:end,conc_ind+1:end,:))],1);
    end
    %variance of mean adjusted contemporaneous concordance index
    [ac, cov, hm, sd, med, ki, ineff, w] = autocovneu_conc(It,10);
    cv = conc_var(cov,size(cov,1)-1,conc_ind,w);

    t=['with series ' char(ctry(conc_ind)) ' \n'];
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
    ind_exconc=[conc_ind+1:size(It,2)];
    for i=1:length(ind_exconc)
        t1=[];
        t2=[];
        t1=[char(ctry(ind_exconc(i))) '&' num2str(conc_mat(1,i,1)./sqrt(cv(i)),'%5.2f') ];
        t2=['unadj. &  ' ];
        for j=1:maxlead+1;
            t1=[ t1 '& ' num2str(conc_mat(j,i,1),'%5.2f') ];
            t2=[ t2 '& ' num2str(conc_mat(j,i,2),'%5.2f') ];
        end
        %         t1 = [t1 '\\\\ \n'];
        fprintf(fid,[t1 '\\\\ \n']);
        fprintf(fid,[t2 '\\\\ \n']);
    end
    t='\\hline \\hline \\end{tabular}\\end{center}\\end{table}\n \n';
    fprintf(fid,t);
end


% fprintf(fid,'P-Value GDP \n');
% % t=num2str([P(3,:);P(end,:)],'%8.4f %8.4f\n');
% t = [P(1,:)];
% fprintf(fid,'%6.4f\n',t);

% figure(gcf+1)
% sig=sign(conc_mat(find(abs(conc_mat)==kron(ones(maxlead+1,1),max(abs(conc_mat),[],1)))));
% scatter(conc_mat(1,:)',sig.*max(abs(conc_mat),[],1)');
% xlabel('concordance index (contemporaneous)')
% ylabel('maximum concordance index')


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

