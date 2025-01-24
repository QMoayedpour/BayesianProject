% fuehrt zu random-permutation sampler im nachhinein eine entsprechende
% Sortierung/Permutation durch

% fuer alphamc, etamc, etapostmc, anmc, ancholmc,(postmc fehlt)

%eta_restr=0;
%group=1;
M=size(alphamc,1);
% mn=size(anmc,1);
etaMSmc=1;
nst=max(1,size(etaMSmc,1));
dd=size(Q,2);
dMS=size(ZMS,2);
MSmc=M-size(Smc,1);
indexMIX=[1:dd*K];
indexFIX=[K*dd+1:d+(K-1)*dd];
indexMS=[(d+(K-1)*dd+1):(d+(K-1)*dd+dMS*K)];
%neu anfang
ZG=Z;
% search_colMS
for i = 1:M
    
    %alphai = reshape(alphamc(i,indexMIX),dd,K);
    alphai=reshape(alphamc(i,[indexMIX indexMS]),dd,K,nst); 
    if nst>1
        etai=etaMSmc(:,:,1,i);
    end
    
    alpha=alphamc(i,:)';
    if dMS>0
        for k=1:K;
            indexMSsort=[1 1 1 1]; %as many elements as groups
            if (alpha(indexMS((k-1)*dMS+indexMSsort(k)))<0)*(1-eta_restr)+(etai(1,1)>etai(2,2))*eta_restr
                %                 if etai(1,1)>etai(2,2) 
                %                             if alphai(indexMSsort(k),k,1)<alphai(indexMSsort(k),k,1)-alphai(indexMSsort(k),k,2)
                % match group specific and fixed parameters 
                iiMIX=index_perm2(index_perm2<=dd);
                %iiFIX=index_perm2(index_perm2>dd);
                alphai(iiMIX,k,1)=alphai(iiMIX,k,1)-alphai(iiMIX,k,2);
                %alpha(indexFIX(iiFIX))=alpha(indexFIX(iiFIX))-alpha(indexMS(iiFIX)); 
                alphai(iiMIX,k,2)=-alphai(iiMIX,k,2);
                if i>MSmc
                    if size(IMSmc,1)>1
                        IMSmc(k,:,i-MSmc) = 1-IMSmc(k,:,i-MSmc);
                    elseif k==1
                        IMSmc(k,:,i-MSmc) = 1-IMSmc(k,:,i-MSmc);
                    end
                end   
                if size(etaMSmc,3)>1
                    etaMSmc(:,:,k,i)=etaMSmc([2 1],[2 1],k,i);
                elseif k==1
                    etaMSmc(:,:,k,i)=etaMSmc([2 1],[2 1],k,i); 
                end
            end
        end
        alphamc(i,indexMIX)=reshape(alphai(:,:,1),dd*K,1)';
        alphamc(i,indexMS)=reshape(alphai(:,:,2),dd*K,1)';
        
    end
    
    % Group specific parameters
    alpha1 = alphai;
    alpha1(dd+1,:,:)=permute(kron(ones(nst,1),[1:K]),[3 2 1]);
    if group
        if  or(group==1, group==2)
            if group==1
                [alph,ssi] = sort(alpha1(2,:,1)');
                alpha1 = alpha1(:,ssi,:);
                %             [alph,ssi] = sort(alpha1(2,1:2,1)');
                %             alpha1(:,1:2,:) = alpha1(:,ssi,:);
            elseif group==2
                [alph,ssi] = sort(alpha1(1,:,2)');
                alpha1 = alpha1(:,ssi,:);
            end
            si = alpha1(dd+1,:,1);
            
        else
            if group==3
                
                if K==2
                    if Qmc(i,1,1)>Qmc(i,1,2)
                        si=[2 1];
                    else
                        si=[1 2];
                    end   
                else
                    'K >2 einbauen'
                end   
            end
        end
        %     [alph,ssi] = sort([alpha1(1,2:3,1)-alpha1(1,2:3,2)]');
        %     alpha1(:,2:3,:)=alpha1(:,ssi+1,:);
        %       [alph1,ssi1] = sort(squeeze(alpha1(1,1:2)));
        %       alpha1(:,1:2) = alpha1(:,ssi1);
        
        % sort such that first group is the largest one
        %[se,ssi]=sort(etamc(i,:));
        %ssi=ssi(end:-1:1);
        %alpha1 = alpha1(:,ssi);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    alphamc(i,indexMIX) = reshape(alphai(:,si,1),dd*K,1)';
    if nst>1
        alphamc(i,indexMS) = reshape(alphai(:,si,2),dd*K,1)';
    end
    
    etamc(i,:) = etamc(i,si);
    if and(nst>1,size(etaMSmc,3)>1);
        etaMSmc(:,:,:,i)=etaMSmc(:,:,si,i);
    end
    Qmc(i,:,:)=Qmc(i,:,si);
    
    if i>MSmc
        [si1,ssort]=sort(si);
        Smc(i-MSmc,:) = ssort(Smc(i-MSmc,:));
        if nst>1
            if size(IMSmc,1)>1;
                IMSmc(:,:,i-MSmc) = IMSmc(si,:,i-MSmc);
            end
        end
    end      
    
    
end      
am_MIX=reshape(mean(alphamc(:,indexMIX)),dd,K);
am_FIX=mean(alphamc(:,indexFIX));
t_MIX=reshape(mean(alphamc(:,indexMIX)),dd,K)./reshape(std(alphamc(:,indexMIX)),dd,K);
t_FIX=mean(alphamc(:,indexFIX))./std(alphamc(:,indexFIX));
sd_MIX=reshape(std(alphamc(:,indexMIX)),dd,K);
sd_FIX=std(alphamc(:,indexFIX));
if nst>1
    t_MS=reshape(mean(alphamc(:,indexMS))./std(alphamc(:,indexMS)),dMS,K);
    am_MS=reshape(mean(alphamc(:,indexMS)),dMS,K);
    sd_MS=reshape(std(alphamc(:,indexMS)),dMS,K);
    
    am_REC=reshape([mean([reshape(alphamc(:,indexMIX),M,dd,K)-reshape(alphamc(:,indexMS),M,dd,K)],1)],dd,K);
    t_REC=reshape([(mean([reshape(alphamc(:,indexMIX),M,dd,K)-reshape(alphamc(:,indexMS),M,dd,K)],1))./(std([reshape(alphamc(:,indexMIX),M,dd,K)-reshape(alphamc(:,indexMS),M,dd,K)],[],1))],dd,K);
    sd_REC=reshape(std([reshape(alphamc(:,indexMIX),M,dd,K)-reshape(alphamc(:,indexMS),M,dd,K)],[],1),dd,K);
    [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc alphamc(:,indexMIX)-alphamc(:,indexMS)],10);
    ki_MIX=reshape(ki(indexMIX,:),dd,K,2);
    ki_REC=reshape(ki(indexMS(end)+1:end,:),dd,K,2);
    ki_MS=reshape(ki(indexMS,:),dd,K,2);
end

if nst==1
    [ac, hm, sd, med, ki, ineff] = autocovneu([alphamc],10);
    ki_MIX=reshape(ki(indexMIX,:),dd,K,2);
end


if nst>1;
    etadiag=[];
    for j=1:size(etaMSmc,1)
        etadiag=[etadiag permute(etaMSmc(j,j,:,:),[4 3 1 2])];
    end
    [ac, hm, sd, med, ki2, ineff] = autocovneu(etadiag,10);
    ki_eta=reshape(ki2(1:size(etaMSmc,3),:),1,size(etaMSmc,3),2);
    for j=2:nst
        ki_eta=[ki_eta; reshape(ki2((j-1)*size(etaMSmc,3)+1:j*size(etaMSmc,3),:),1,size(etaMSmc,3),2)];
    end
end
% am_MIXSUM=squeeze(mean(sum(reshape(alphamc(:,indexMIX),M,dd,K),2),1))
% sd_MIXSUM=squeeze(std(sum(reshape(alphamc(:,indexMIX),M,dd,K),2)))
% t_MIXSUM=squeeze(mean(sum(reshape(alphamc(:,indexMIX),M,dd,K),2),1))./squeeze(std(sum(reshape(alphamc(:,indexMIX),M,dd,K),2)))
% am_MSSUM=squeeze(mean(sum([reshape(alphamc(:,indexMIX),M,dd,K)-alphamc(:,indexMS,ones(1,K))],2),1))
% sd_MSSUM=squeeze(std(sum([reshape(alphamc(:,indexMIX),M,dd,K)-alphamc(:,indexMS,ones(1,K))],2)))
% t_MSSUM=squeeze(mean(sum([reshape(alphamc(:,indexMIX),M,dd,K)-alphamc(:,indexMS,ones(1,K))],2),1))./squeeze(std(sum([reshape(alphamc(:,indexMIX),M,dd,K)-alphamc(:,indexMS,ones(1,K))],2)))

cov_MIXMS=cov([alphamc(:,indexMIX) alphamc(:,indexMS)]);
if nst>1
    m_eta=mean(etaMSmc,4);
    s_eta=sort(etaMSmc,4);
end