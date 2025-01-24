
% [data ser]=xlsread('DAT_roh_Q0410.xls'); %
[data ser]=xlsread('DAT_roh_Q0604.xls'); % data partially available till 2006Q1
% the first line, 1987Q4 contains an indicator which relates to a specific class of series (see printout)
% 1   GDP and its main components
% 2   WIFO business cycle indicator
% 3   business surveys
% 4   prices (HICP)
% 5   CPI
% 6   wages
% 7   wholesale prices
% 8   foreign trade
% 9   labor market
% 10  IP
% 11  financial variables

% Zlogit=load('corr_0604.asc');

% ser_excl={'SZR';'YED';'PCD';'ITD';'GCD';'MTD';'XTD'}; %string of series to exclude
ser_excl={'SZR';'YED';'PCD';'ITD';'GCD';'MTD';'XTD';'ECOSEN';'INDSEN'}; %string of series to exclude
ser_diff={'SZR';'QTAUF';'QTEXPA';'QTLAG';'QTPR';'QTPRO';'QTBAUF';'QTBPR';'QTBBGL';'QTBAGL';'INDSEN';'KTPROL';'KTAUF';'KTAUSL';...
            'KTLAG';'KTPRON';'KTVPN';'BAUVPN';'EINDSE';'EBAUSE';'EHANSE';'EKONSE';'ALQN';'ALQNSA';'STI';'SEKMRE';'YIELD'};
% ser_cont={'YER';'PCR';'ITR';'MTR';'XTR'}; %series falling into the contemporaneous group
ser_cont={'YER';'PCR';'ITR';'MTR';'XTR'}; %series falling into the contemporaneous group
ser_lead={'QTAUF';'KTPROL';'KTAUF';'EINDSE'}; %series falling into the leading group
% ser_incl=[ser_cont;ser_lead;{'QTBAUF';'KTPRON';'EXP7';'DAX';'YIELD';'DEBT';'GHPIVBG';'VPIG86'}];

ser=ser(1,2:end); % string array of data series
% ser_excl=ser([4 7:end]); %string of series to exclude
ser_sel = ser(find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1))); % string array of included series
ser_ind=find(1-sum(strcmp(ser(ones(size(ser_excl,1),1),:),ser_excl(:,ones(1,size(ser,2)))),1)); %choose which countries to include, index array
% ser_sel = ser(find(sum(strcmp(ser(ones(size(ser_incl,1),1),:),ser_incl(:,ones(1,size(ser,2)))),1))); % string array of included series
% ser_ind=find(sum(strcmp(ser(ones(size(ser_incl,1),1),:),ser_incl(:,ones(1,size(ser,2)))),1)); %choose which countries to include, index array
ind_cont=find(sum(strcmp(ser_sel(ones(size(ser_cont,1),1),:),ser_cont(:,ones(1,size(ser_sel,2)))),1)); %index of included series defining the contemporaneous group
ind_lead=find(sum(strcmp(ser_sel(ones(size(ser_lead,1),1),:),ser_lead(:,ones(1,size(ser_sel,2)))),1)); %index of included series defining the leading group

cat_dum=data(1,:);
cat_dum=cat_dum(ser_ind);
data=data(2:end,:);
% data=log(data);
% y_gr=(data(2:end,:)-data(1:end-1,:))*100;
% y_gr=y_gr- kron(ones(size(y_gr,1),1),mean(y_gr,1));
% cctries=cellstr(ctries_sel);
ind_diff=find(sum(strcmp(ser(ones(size(ser_diff,1),1),:),ser_diff(:,ones(1,size(ser,2)))),1));
ind_logdiff=find(1-sum(strcmp(ser(ones(size(ser_diff,1),1),:),ser_diff(:,ones(1,size(ser,2)))),1));
ind_oel=find(strcmp(ser,'OEL'));
ind_mult100=ind_logdiff([1:find(ind_logdiff==ind_oel)-1 find(ind_logdiff==ind_oel)+1:end]);

data(:,ind_logdiff)=log(data(:,ind_logdiff));
data=data(2:end,:)-data(1:end-1,:);
data(:,ind_mult100)=data(:,ind_mult100)*100;



nobs=size(data,1);
% cend = 2004.75; %with DAT_roh_Q0410.xls
cend = 2006.25; %with DAT_roh_Q0412.xls
cal = [sort(cend-[1:nobs]/4)]';
% cal_beg=cal(1);cal_end=cal(end);

% lrmean_smpl=find((cal>=cal_beg).*(cal<cal_end+1)); %choose sample period to compute the long-run mean of the series
obs_per=find((cal>=cal_beg).*(cal<cal_end+1)); %choose sample period, contains the respective rows of original matrix
% lrmean_smpl=find((cal>=1980).*(cal<2003)); %choose sample period to compute the long-run mean of the series
% obs_per=find((cal>=1980).*(cal<2003)); %choose sample period, contains the respective rows of original matrix
cal=cal(obs_per);

nbank=length(ser_ind);
S=zeros(1,nbank);


% dlo=(log(data(2:end,:))-log(data(1:end-1,:)))*100; %endogenous variable, GDP quarterly growth rate
dlo=data(:,ser_ind);
% Zlogit=Zlogit(ser_ind,:);

dlo=dlo(obs_per(1:end),:); %ev. end-1 if data is differenced on line 51
dlo_unadj=dlo;
nobs=size(dlo,1);
% cal=cal(2:end); %ev. uncomment if data is differenced on line 51
% dlo_obs=dlo;

%set outliers
outl=zeros(nobs,nbank);
ser_outl={'PCR';'GCR';'MTR'};
J=[find(sum(strcmp(ser_sel(ones(size(ser_outl,1),1),:),ser_outl(:,ones(1,size(ser_sel,2)))),1))];
for j=J;
    if j==J(1)
        %         r=[find(cal==1996.25),find(cal==1996.50)];
        r=[find(cal==1993)];

        outl(r,j)=1;
    end
    if j==J(2)
        %         r=[find(cal==1994.75),find(cal==1995),find(cal==1997)];
        r=[find(cal==1996)];
        outl(r,j)=1;
    end
    if j==J(3)
        %         r=[find(cal==1993.0),find(cal==1993.25)];
        r=[find(cal==1993)];
        outl(r,j)=1;
    end
end


nan_ret=find(1-((mean(dlo,1)>-Inf).*(mean(dlo,1)<Inf))); %find series with NaN  
n_ret=zeros(nobs,nbank);
for j=nan_ret;
   n_ret([find(1-((dlo(:,j)>-Inf).*(dlo(:,j)<Inf)))],j)=1; %construct matrix of missing values for NaN
end

dlo_out=(1-(1-n_ret).*(1-outl))==1;

%compute mean and std without outliers
for j=1:nbank
    lrmean(j)=mean(dlo(dlo_out(:,j)==0,j));
    lrstd(j)=std(dlo(dlo_out(:,j)==0,j));
%    dlo([find(dlo_out(:,j))],j)=mean(dlo(dlo_out(:,j)==0,j));
end

% lrmean=mean(dlo); %ev. end-1 if data is differenced on line 51
% lrstd=std(dlo,0,1);

if ~unit_spec_var
    dlo=(dlo-kron(ones(nobs,1),lrmean))./kron(ones(nobs,1),lrstd);
elseif unit_spec_var
    dlo=(dlo-kron(ones(nobs,1),lrmean));
%    dlo=dlo-kron(ones(nobs,1),lrmean);
end    
dlo_obs=dlo;

if sg_change
    % [sg] = max_conc_sign(dlo,ser_sel,{'YER'},4); %change sign for countercyclical variables according to the sign of largest concordance
    
    [sg] = neg_corr_sign(dlo,ser_sel,{'YER'}); %change sign for series significantly negatively correlated with GDP
    dlo=dlo.*sg(ones(nobs,1),:);
    dlo_obs=dlo_obs.*sg(ones(nobs,1),:);
    dlo_unadj=dlo_unadj.*sg(ones(nobs,1),:);
    [sg] = neg_corr_sign(dlo,ser_sel,{'KTAUF'}); %change sign for series significantly negatively correlated with KTAUF
    dlo=dlo.*sg(ones(nobs,1),:);
    dlo_obs=dlo_obs.*sg(ones(nobs,1),:);
    dlo_unadj=dlo_unadj.*sg(ones(nobs,1),:);

end

if S_logit
    corr_c={'YER'}; %correlation with GDP
    c_ind=find(sum(strcmp(ser_sel(ones(size(corr_c,1),1),:),corr_c(:,ones(1,size(ser_sel,2)))),1));
    corr_l={'KTAUF'}; %correlation with orders
    l_ind=find(sum(strcmp(ser_sel(ones(size(corr_l,1),1),:),corr_l(:,ones(1,size(ser_sel,2)))),1));
    [R P]=corrcoef(dlo);
    Zlogit=[R(c_ind,:)' R(l_ind,:)'];
end



% Zbasis=zeros(nobs,[],nbank);
% Zbasis(:,1,:)=reshape(lenddata_b(:,9),nobs,nbank)*100; % zinssatz dir
% Zbasis(:,2,:)=reshape(lenddata_b(:,3),nobs,nbank); %Grösse si

% Zdum=zeros(nobs,4,nbank);%dummies 1-3: saisonals 4: structural break in der liquidität ab 4.quartal 1995 
% for i1=1:4;
%    Zdum(:,i1,:)=reshape(lenddata_b(:,4+i1),nobs,nbank); 
% end;clear i1;   

% Zwirt=zeros(nobs,2,nbank); % wirtschaftsdaten 1: dp86 (inflation) 2: dyr (wachstum bip)
% for i1=1:2;
%    Zwirt(:,i1,:)=reshape(lenddata_b(:,9+i1),nobs,nbank)*100; 
% end;clear i1;   
% dyr=lenddata_b(1:nobs,9+2)*100;
c=ones(size(dlo,1),1,size(dlo,2)); % constante 


clear outl n_ret 
