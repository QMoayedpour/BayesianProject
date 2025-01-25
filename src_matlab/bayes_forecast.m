function [bfore_mat,bforeIMS_mat]=bayes_forecast(H,Hf,cal,M,vars,fully_b,stdep,y,Z,alphamc,Qmc,sepsmc,lambmc,etaMSmc,IMSmc,ISmc,pool,unit_spec_var,K,dd,d,dMS,isp_lag_dlo);

%T=H;
N=size(ISmc,2);
MC=size(ISmc,3);
nvars=length(vars);
bfore_mat=zeros(nvars,Hf,M,H+1); %H out-of-sample forecasts for evaluation, 1 block with out-of-sample
% forecast at the actual end
bforeIMS_mat=(zeros(1+(K>1),Hf,M,H+1)==1);

if fully_b
    u = simuni(MC,M);
    alpha=alphamc(u,:);
    seps=sepsmc(u,:);
    IMS=IMSmc(:,:,u);
    IS=ISmc(:,:,u);
    etaMS=etaMSmc(:,:,u);
    if unit_spec_var
        lamb=lambmc(u,:);
        seps=seps(:,ones(1,N))./lamb;
    end
    if any(Qmc)
        Q=Qmc(u,:,:);
    else Q=0;
    end
elseif ~fully_b
    alpha=kron(ones(M,1),mean(alphamc,1));
    seps=kron(ones(M,1),mean(sepsmc,1));
    IMS=mean(IMSmc,3)>0.5;
    IMS=IMS(:,:,ones(M,1));
    IS=mean(ISmc,3)==kron(ones(K,1), max(mean(ISmc,3)));
    IS=IS(:,:,ones(M,1));
    etaMS=mean(etaMSmc,4);
    etaMS=etaMS(:,:,ones(M,1));
    if unit_spec_var
        lamb=kron(ones(M,1),mean(lambmc,1));
        seps=seps(:,ones(1,N))./lamb;
    end
    if any(Qmc)
        Q=mean(Qmc,1);
        Q=Q(ones(1,M),:,:);
    else Q=0;
    end
end




%%  define current and future  Z
% run_model_lagzins_prognose;

% % define banks to be consider for forecasting
% fc_beg=1999; % forecast interval begin
% fc_end=2001.75; % forecast interval end
% indexfc = define_banks(fc_beg,fc_end); %adjust procedure, in particular to enforce some banks to be considered despite outliers.
% index_forecast =  find(sum(kron(ones(length(indexfc),1),indexb)'==kron(ones(length(indexb),1),indexfc),2))';
% %adjust index_forecast to the dimension of estimation sample
%
% cwd = pwd;
% cd(tempdir);
% pack
% cd(cwd)


% Zh_m=zeros(Hf,size(Z,2),M);
% Zh_m=Z;
% %%TEST: future y  exogen: put comment on next line
% Zh_m(T+1:T+h,isp_lag_dlo,:)=0; % the remaining variables are assumed to be known
%
% Ma=size(alphamc,1);
% MS=size(Smc,1);
% yh_m=zeros(h,N,MS);
% IMSpro=zeros(h,MS);

% d = size(Z,2); % davon sind die ersten dd Spalten variabel
% dd=lag_dir;dMS=dd;

indexMIX=1:(K)*dd;
indexFIX=[K*dd+1:d+(K-1)*dd];
indexMS=[(d+(K-1)*dd+1):(d+(K-1)*dd+dMS*(K))];
SG=kron(ones(1,N),[1:K]');

for h=1:H+1
    
    for m=1:M
        if h<H+1
            Zf = permute(Z(h,:,:),[2 3 1]);
        elseif h==H+1
            Zf = [permute(Z(end,1,:),[2 3 1]); y(end,:); permute(Z(end,2:end-1,:),[2 3 1])]; %adjust in case of exogenous variables.
        end
        IMSf = IMS(:,h,m);
        yf = zeros(1,N);
        err = 0;

        S = SG(IS(:,:,m));
        betaG=zeros(dd,N);
        betaR=zeros(dMS,N);

        for   k=1:1+(K>1)
            betaG(:,S==k)=alpha(m*ones(sum(S==k),1),indexMIX((k-1)*dd+[1:dd]))';
            if dMS>0
                betaR(:,S==k)=alpha(m*ones(sum(S==k),1),indexMS((k-1)*dMS+[1:dMS]))';
            end
        end
        beta=[betaG;alpha(m*ones(N,1),indexFIX)'];

        for hf=1:Hf;

            %         Zh_m(t,isp_lag_dlo(2:end),:)=Zh_m(t-1,isp_lag_dlo(1:end-1),:);
            %         if t>T+1
            %             %%TEST: future y  exogen: put comment on next line
            %             Zh_m(t,isp_lag_dlo(1),:)=yh_m(t-T-1,:,m);
            if dMS>0
                IMSf = forecast_MS(IMSf(1:1+(K>1),:),etaMS(:,:,m));
                if K>1
                IMSf = [IMSf;1];
                end
                
            end

            %    yh_m(t-T,:,m)=squeeze(Zh_m(t,1,:))'.*alphamc(Ma-MS+m,Smc(m,:))+squeeze(alphamc(Ma-MS+m,end-1:end))*squeeze(Zh_m(t,[2 3],:));
%            if fully_b
                err= seps(m,:).^.5;
                if any(Q);
                    for k=1:1+(K>1)
                        kind=find(S==k);
                        ksum=sum(S==k);
                        ZQ=[Zf(:,kind)' Zf(1:dMS,kind)'.*(kron(IMSf(k),ones(ksum,dMS))-1)];
                        Qi=qinmatr(Q(m,:,k)');
                        err(kind)= (diag(ZQ*Qi*ZQ').^.5)' + err(kind);
                    end
                end
                err=err.*randn(1,N);
%            end
            betaall=beta;
            if dMS>0
                betaall(1:dMS,:)=betaall(1:dMS,:)+ (kron(ones(dMS,1),reshape(IMSf(S),1,N))-1).*betaR;
            end
            yf = err +  sum(Zf.*betaall,1);
            bfore_mat(:,hf,m,h)= yf(1,vars);
            bforeIMS_mat(:,hf,m,h)=IMSf(1:1+(K>1));
            Zf=[Zf(1,:); yf; Zf(2:end-1,:)] ; %adjust in case of exogenous variables

        end

    end
end

% ytrue=y(T+1:T+h,:);
% errall= (yh_m(:,index_forecast,:)-ytrue(:,index_forecast,ones(1,MS)));
%
% w=lambmc(Ma-MS+[1:MS],index_forecast).^.5;
% errall_w=errall.*permute(w(:,:,ones(h,1)),[3 2 1]);
%
% w_si=kron(ones(MS,1),exp(si(cal==1998.75,index_forecast))./sum(exp(si(cal==1999.75,index_forecast))));
% errall_wsi=errall.*permute(w_si(:,:,ones(h,1)),[3 2 1]);
%
%
% biasall(imodel,:,1)=mean(mean(errall,2),3)';
% mseall(imodel,:,1)=mean(mean(errall.^2,2),3)';
% biasall(imodel,:,2)=mean(mean(errall_w,2),3)';
% mseall(imodel,:,2)=mean(mean(errall_w.^2,2),3)';
% biasall(imodel,:,3)=mean(mean(errall_wsi,2),3)'
% mseall(imodel,:,3)=mean(mean(errall_wsi.^2,2),3)'




