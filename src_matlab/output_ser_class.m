fid=fopen("fileout2.out",'a');
fprintf(fid,'\n \n');
K=size(etamc,2);
t=['\\begin{table}\\caption{Series classification}\\label{tab:serclass}\\begin{center}\\begin{tabular}{|' ... 
        '*{' int2str(K) '}{l}|} \\hline \\hline \n'];
fprintf(fid,t);

rho=mean(rhomc)';

t1=['$S_i=1$ '];
t2=['$\\rho$=' num2str(rho(1),'%5.2f') ' '];
for j=2:K;
    t1=[t1 '& $S_i=' int2str(j) '$'];
    t2=[t2 '& $\\rho$=' num2str(rho(j),'%5.2f')];
end
t1=[t1 ' \\\\ \n'];
t2=[t2 ' \\\\ \\hline \n'];

fprintf(fid,t1);
fprintf(fid,t2);

n_grp=sum(prob,2);
max_ngrp=max(n_grp);
emp=' ';
amper='&';
amper=amper(ones(max_ngrp,1));
endline=' \\\\ ';
endline=endline(ones(max_ngrp,1),:);
t= [strvcat(ser_sel(prob(1,:)==1)) ; emp(ones(max_ngrp-n_grp(1),1),ones(size(strvcat(ser_sel(prob(1,:)==1)),2),1))];
for j=2:K;
    t=[t amper [strvcat(ser_sel(prob(j,:)==1)) ; emp(ones(max_ngrp-n_grp(j),1),ones(size(strvcat(ser_sel(prob(j,:)==1)),2),1))]];    
end
t=[t endline];
for j=1:max_ngrp;
    txt=[t(j,:) '\n'];
   fprintf(fid,[t(j,:) '\n']);
end
t=['\\\\ \\hline \\hline \\end{tabular}\\end{center}\\end{table} \n'];
fprintf(fid,t);
% 
% % t2=[];
% % t3=[];
% for i=1:N;
%     t1=[];
%     t1=[char(ctry(i)) ];
%        for j=1:size(prop_mat,1)-2;
%            t1=[ t1 '& ' num2str(prop_mat(j,i),'%5.2f') ];
%        end
%        if P(1,i)<0.05
%            t1=[ t1 '& {\\bf ' num2str(prop_mat(j+1,i),'%5.2f') '}' ];
%        else
%            t1=[ t1 '& ' num2str(prop_mat(j+1,i),'%5.2f') ];
%        end
%        t1=[ t1 '& ' num2str(prop_mat(end,i),'%5.2f') ];
%        t1 = [t1 '\\\\ \n'];
%        fprintf(fid,t1);
% end
% t='\\hline \\hline \\end{tabular}\\end{center}\\end{table}\n';
% fprintf(fid,t);
% fprintf(fid,'P-Value GDP \n');
% % t=num2str([P(3,:);P(end,:)],'%8.4f %8.4f\n');
% t = [P(1,:)];
% fprintf(fid,'%6.4f\n',t);
% 
% % fprintf(fid,'80-02   80-87  85-92  90-97  95-02\n');
% % fprintf(fid,'%6.2f  %6.2f %6.2f %6.2f %6.2f\n',mean_mat);
% 
% % t=num2str(P(end,:),'%8.4f\n');
% % fprintf(fid,'P-Value Corr with US \n');
% 
% 
% % figure(gcf+1)
% % scatter(R8092(3,[1:2 4:end-1]),R8092(end,[1:2 4:end-1]),'o','filled');
% % % scatter(R8092(3,:),R8092(end,:));
% % title('correlations of GDP growth 1980-1992')
% % xlabel('correlation with German GDP')
% % ylabel('correlation with US GDP')
% % hold on;
% % for j=[1 2 4:N-1];
% % % for j=1:N;
% %     text(R8092(3,j)-0.01,R8092(end,j),ctries_sel(j,:),'VerticalAlignment','middle','HorizontalAlignment','right','FontSize',12)
% % end
% 
% % figure(gcf+1)
% % scatter(R9002(3,[1:2 4:end-1]),R9002(end,[1:2 4:end-1]),'o','filled');
% % % scatter(R8092(3,:),R8092(end,:));
% % title('correlations of GDP growth 1990-2002')
% % xlabel('correlation with German GDP')
% % ylabel('correlation with US GDP')
% % hold on;
% % for j=[1 2 4:N-1];
% % % for j=1:N;
% %     text(R9002(3,j)-0.01,R9002(end,j),ctries_sel(j,:),'VerticalAlignment','middle','HorizontalAlignment','right','FontSize',12)
% % end
% 
% % 
% 
% % [corr_sor si]=sort(diag(Rlag(N+1:end,1:N)));
% % lab_sor=ctries_sel(si,:);
% % figure(gcf+1)
% % bar(1:N,corr_sor);
% % % scatter(R8092(3,:),R8092(end,:));
% % title('correlation of y_{it} with y_{i,t-1}')
% % ylabel('correlation')
% % hold on;
% % for j=[1:N];
% % % for j=1:N;
% %     text(j,(corr_sor(1)<0)*(corr_sor(1)-0.07),lab_sor(j,:),...
% %                'VerticalAlignment','middle','HorizontalAlignment','right','Rotation',90,'FontSize',12)
% % end
% % 
% % t=[ctries_sel num2str(diag(Rlag(N+1:end,1:N)),'%8.4f')  num2str(diag(Plag(N+1:end,1:N)),'%8.4f')]
st=fclose(fid);
% % 
% % figure(gcf+1)
% % for j=[1:N];
% % scatter(y_gr(1:end-1,j),y_gr(2:end,j));
% % % scatter(R8092(3,:),R8092(end,:));
% % hold on;
% % end
% % title('y_{t-1} against y_t')
% % ylabel('y_{t}')
% % xlabel('y_{t-1}')
