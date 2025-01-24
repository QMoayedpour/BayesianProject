%ind_enc=[2 3]; %first contemporaneous, second leading
%file=eingeben;
t_val=2; %flag for significance (0 for s.e.,1 for t-value,2 for shortest 95% conf.int.), 
euro=1; %flag to convert to euro

nlag=size(am_MIX,1);
ngrp=size(am_MIX,2);

permut_MS_ALL;
fid=fopen('test.out','w');
t=['\\begin{table}\\caption{ }\\label{tab:1}\\begin{center}\\begin{tabular}{l|' ... 
        '*{' int2str(ngrp) '}{c}|*{' int2str(ngrp) '}{c}|*{' int2str(ngrp) '}{c}|} \\hline \\hline \n'];
fprintf(fid,t);

t=['coeff.& \\multicolumn{' int2str(ngrp) '}{l|}{$I_t=1$}& \\multicolumn{' int2str(ngrp*2) '}{l|}{$I_t=0$}  \\\\ \n'];
fprintf(fid,t);
t1=[];
t2=[];
t3=[];
for i=1:ngrp;
    t1=[t1 '&$S_i=' int2str(i) '$'];
    t2=[t2 '&$S_i=' int2str(i) '$'];
    t3=[t3 '&$\\mu^R_{' int2str(i) '\\cdot}$' ];
end
t=[t1 t2 t3 ' \\\\ \\hline \n'];

fprintf(fid,t);

for j=1:nlag;
    t11=['coeff$_' int2str(j) '$&' num2str(am_MIX(j,1),'%5.2f') '&' ];
    if j<=dMS;
        t12=[num2str(am_REC(j,1),'%5.2f') '&'  ];
        t13=[num2str(am_MS(j,1),'%5.2f') ];
    end
    if t_val==1
        t11b=[ '& (' num2str(t_MIX(j,1),'%5.2f') ')&'];
        if j<=dMS;
            t12b=[ '(' num2str(t_REC(j,1),'%5.2f') ')&'];
            t13b=[ '(' num2str(t_MS(j,1),'%5.2f') ')'];
        end
    elseif t_val==0
        t11b=[ '& (' num2str(sd_MIX(j,1),'%5.2f') ')&'];
        if j<=dMS;
            t12b=[ '(' num2str(sd_REC(j,1),'%5.2f') ')&'];
            t13b=[ '(' num2str(sd_MS(j,1),'%5.2f') ')'];
        end
    elseif t_val==2
        t11b=[ '& (' num2str(ki_MIX(j,1,1),'%5.2f') ' ' num2str(ki_MIX(j,1,2),'%5.2f') ')&'];
        if j<=dMS
            t12b=[ '(' num2str(ki_REC(j,1,1),'%5.2f') ' ' num2str(ki_REC(j,1,2),'%5.2f') ')&'];
            t13b=[ '(' num2str(ki_MS(j,1,1),'%5.2f') ' ' num2str(ki_MS(j,1,2),'%5.2f') ')'];
        end
    end
    for i=2:ngrp;
        t11=[t11 num2str(am_MIX(j,i),'%5.2f') '&' ] ;
        if j<=dMS
            t12=[t12 num2str(am_REC(j,i),'%5.2f') '&'  ];
            t13=[t13 '&' num2str(am_MS(j,i),'%5.2f') ];
        end
        if t_val==1
            t11b=[t11b '(' num2str(t_MIX(j,i),'%5.2f') ')&'];
            if j<=dMS
                t12b=[t12b '(' num2str(t_REC(j,i),'%5.2f') ')&'];   
                t13b=[t13b '&(' num2str(t_MS(j,i),'%5.2f') ')'];
            end
        elseif t_val==0
            t11b=[t11b '(' num2str(sd_MIX(j,i),'%5.2f') ')&'];
            if j<=dMS
                t12b=[t12b '(' num2str(sd_REC(j,i),'%5.2f') ')&'];  
                t13b=[t13b '&(' num2str(sd_MS(j,i),'%5.2f') ')'];
            end
        elseif t_val==2
            t11b=[t11b '(' num2str(ki_MIX(j,i,1),'%5.2f') ' ' num2str(ki_MIX(j,i,2),'%5.2f') ')&'];
            if j<=dMS
                t12b=[t12b '(' num2str(ki_REC(j,i,1),'%5.2f') ' ' num2str(ki_REC(j,i,2),'%5.2f') ')&'];
                t13b=[t13b '&(' num2str(ki_MS(j,i,1),'%5.2f') ' ' num2str(ki_MS(j,i,2),'%5.2f') ')'];
            end
        end
    end
    t1=[t11] ;
    t2=[t11b];
    if j<=dMS
        t1=[t1 t12 t13];
        t2 =[t2 t12b t13b];
    elseif j>dMS
        t1=[t1 '\\multicolumn{' int2str(ngrp*(2*(nst>1))) '}{l|}{\\ }' ];
        t2 =[t2 '\\multicolumn{' int2str(ngrp*(2*(nst>1))) '}{l|}{\\ }'];
        
    end
    t1=[t1 ' \\\\ \n'];
    t2=[t2 ' \\\\ \n'];
    fprintf(fid,t1);
    fprintf(fid,t2);
end
fprintf(fid,'\\hline \n');
% t11=[];t12=[];t11b=[];t12b=[];
% for i=1:ngrp;
%     t11=[t11 num2str(am_MIXSUM(i),'%5.2f') '&' ] ;
%     t12=[t12 num2str(am_MSSUM(i),'%5.2f') '&'  ];
%     if t_val
%         t11b=[t11b '(' num2str(t_MIXSUM(i),'%5.2f') ')&'];
%         t12b=[t12b '(' num2str(t_MSSUM(i),'%5.2f') ')&'];   
%     elseif ~t_val
%         t11b=[t11b '(' num2str(sd_MIXSUM(i),'%5.2f') ')&'];
%         t12b=[t12b '(' num2str(sd_MSSUM(i),'%5.2f') ')&'];   
%     end
% end
% t1=['sum&' t11 t12 ];
% t2=['&' t11b t12b ];
% t1=[t1 ' \\\\ \n'];
% t2=[t2 ' \\\\ \n'];
% fprintf(fid,t1);
% fprintf(fid,t2);

if K>=1;
    t11=[];%t12=[];
    t11b=[];%t12b=[];
    t13=[];%t14=[];
    for i=1:ngrp;
        t11=[t11 num2str(am_unc(1,i),'%5.2f') '&' ] ;
        %         t12=[t12 num2str(mean(Zbasis(24,3,class(:,i)==1)),'%5.2f') '&'  ];
        t11b=[t11b '(' num2str(ki_unc(1,i,1),'%5.2f') ' ' num2str(ki_unc(1,i,2),'%5.2f') ')&'];
        %         t12b=[t12b num2str(median(Zbasis(24,3,class(:,i)==1)),'%5.2f') '&'];   
        t13=[t13 num2str(sum(prob(i,:),2)','%5.0f') '&'];
        %     t14=[t14 num2str(sum(exp(loans(size(Zbasis,1)+1,class(:,i)==1)))/sum(exp(loans(size(Zbasis,1)+1,any(class(:,1:ngrp),2)))),'%5.3f') '&'];
        
    end
    if nst>1
        for i=1:ngrp;
            t11=[t11 num2str(am_unc_rec(1,i),'%5.2f') '&' ] ;
            %         t12=[t12 num2str(mean(Zbasis(24,3,class(:,i)==1)),'%5.2f') '&'  ];
            t11b=[t11b '(' num2str(ki_unc(1,ngrp+i,1),'%5.2f') ' ' num2str(ki_unc(1,ngrp+i,2),'%5.2f') ')&'];
            %         t12b=[t12b num2str(median(Zbasis(24,3,class(:,i)==1)),'%5.2f') '&'];   
            %         t13=[t13 num2str(sum(prob(i,:),2)','%5.0f') '&'];
            %     t14=[t14 num2str(sum(exp(loans(size(Zbasis,1)+1,class(:,i)==1)))/sum(exp(loans(size(Zbasis,1)+1,any(class(:,1:ngrp),2)))),'%5.3f') '&'];
        end
        
    end
    t11=['unc.mean&' t11 '\\multicolumn{' int2str(ngrp*(1-(nst==1))) '}{l|}{\\ } \\\\ \n'];
    % t12=['av. liq.&' t12 '\\multicolumn{' int2str(ngrp+1) '}{l|}{\\ } \\\\ \n'];
    t11b=['&' t11b '\\multicolumn{' int2str(ngrp*(1-(nst==1))) '}{l|}{\\ } \\\\  \\hline \\hline \n'];
    % t12b=['med. liq&' t12b '\\multicolumn{' int2str(ngrp+1) '}{l|}{\\ } \\\\ \n'];
    t13=['no. units&' t13 '\\multicolumn{' int2str(ngrp*(2*(nst>1))) '}{l|}{\\ } \\\\ \n'];
    % t14=['market loan sh.&' t14 '\\multicolumn{' int2str(ngrp+1) '}{l|}{\\ } \\\\ \n'];
    fprintf(fid,t11);
    fprintf(fid,t11b);
    fprintf(fid,t13);
    % fprintf(fid,t12);
    % fprintf(fid,t12b);
    % fprintf(fid,t14);
end

if nst>1
    t11=[];t11b=[];t11c=[];t12=[];t12b=[];t12c=[];;
    %t13=[];t14=[];
%    for i=1:size(m_eta,1);
        t11=[diag(m_eta)';ki_eta] ;
%     t11=['$\\xi_{jj}$&' t11 '\\multicolumn{' int2str(ngrp*(2*(nst>1))) '}{l|}{\\ } \\\\ \n'];
%     t12=['$\\xi_{11}$&' t12 '\\multicolumn{' int2str(ngrp*(2*(nst>1))) '}{l|}{\\ } \\\\ \n'];
%     t11b=['conf.int.&' t11b '\\multicolumn{' int2str(ngrp*(2*(nst>1))) '}{l|}{\\ } \\\\ \n'];
%     t12b=['conf.int.&' t12b '\\multicolumn{' int2str(ngrp*(2*(nst>1))) '}{l|}{\\ } \\\\ \n'];
%     t11c=['quarters&' t11c '\\multicolumn{' int2str(ngrp*(2*(nst>1))) '}{l|}{\\ } \\\\ \n'];
%     t12c=['quarters&' t12c '\\multicolumn{' int2str(ngrp*(2*(nst>1))) '}{l|}{\\ } \\\\ \n'];
    %     t13=['no. countries&' t13 '\\multicolumn{' int2str(ngrp+1) '}{l|}{\\ } \\\\ \n'];
    % t14=['market loan sh.&' t14 '\\multicolumn{' int2str(ngrp+1) '}{l|}{\\ } \\\\ \n'];
    %     fprintf(fid,t13);
    fprintf(fid,'%5.2f  %5.2f %5.2f \n',t11);
%     fprintf(fid,t11b);
%     fprintf(fid,t11c);
%     fprintf(fid,t12);
%     fprintf(fid,t12b);
%     fprintf(fid,t12c);
    % fprintf(fid,t14);
    
    t='\\hline \\hline \\end{tabular}\\end{center}\\end{table}\n';
    fprintf(fid,t);
end

 fprintf(fid,'\n');
 fprintf(fid,'persistences \n');
 fprintf(fid,'%5.2f\n',mean(1./(1-etadiag)));
% fprintf(fid,'lead out of trough \n');
% t=[num2str(1/(1-etadiag_enc_m(2)),'%5.2f') '\n'];
% fprintf(fid,t);
% fprintf(fid,'lead of peak \n');
% t=[num2str(1/(1-etadiag_enc_m(3)),'%5.2f') '\n'];
% fprintf(fid,t);


st=fclose(fid);
%[ac, hm, sd, ki, ineff] = autocovneu(par,1000);


