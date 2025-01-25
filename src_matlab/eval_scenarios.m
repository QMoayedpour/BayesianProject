%fid=fopen('eval_scen.txt','w');
load scenarios

scen={'i';'ii';'iii';'iv';'v';'vi';'vii';'viii'};

T=length(i_2003_si);
G=length(i_2003_ci);

It_g=cat(3,[ i_2003_ci i_2003_li],[ i_2006_ci i_2006_li])>0.5;
It_g=cat(2,It_g,cat(3,[ ii_2003_ci ii_2003_li],[ii_2006_ci ii_2006_li])>0.5);
It_g=cat(2,It_g,cat(3,[ iii_2003_ci iii_2003_li],[ iii_2006_ci iii_2006_li])>0.5);
It_g=cat(2,It_g,cat(3,[ iv_2003_ci iv_2003_li],[iv_2006_ci iv_2006_li])>0.5);
It_g=cat(2,It_g,cat(3,[ v_2003_ci v_2003_li],[ v_2006_ci v_2006_li])>0.5);
It_g=cat(2,It_g,cat(3,[ vi_2003_ci vi_2003_li],[ vi_2006_ci vi_2006_li])>0.5);
It_g=cat(2,It_g,cat(3,[ vii_2003_ci vii_2003_li],[ vii_2006_ci vii_2006_li])>0.5);
It_g=cat(2,It_g,cat(3,[ viii_2003_ci viii_2003_li],[ viii_2006_ci viii_2006_li])>0.5);


It=cat(3,[i_2003_si ],[i_2006_si(1:T) ])>0.5;
It=cat(2,It,cat(3,[ii_2003_si ],[ii_2006_si(1:T) ])>0.5);
It=cat(2,It,cat(3,[iii_2003_si],[iii_2006_si(1:T)])>0.5);
It=cat(2,It,cat(3,[iv_2003_si ],[iv_2006_si(1:T) ])>0.5);
It=cat(2,It,cat(3,[v_2003_si ],[v_2006_si(1:T)])>0.5);
It=cat(2,It,cat(3,[vi_2003_si],[vi_2006_si(1:T)])>0.5);
It=cat(2,It,cat(3,[vii_2003_si],[vii_2006_si(1:T)])>0.5);
It=cat(2,It,cat(3,[viii_2003_si],[viii_2006_si(1:T)])>0.5);
N=size(It_g,2)+size(It,2);
maxlead=0;
sIt=size(It,2);
sItg=size(It_g,2);


%concordance matrix : (0,..maxlead) x N x (cyclical, countercyclical)

conc_mat=zeros(maxlead+1,N);

for l=0:maxlead
    %conc_mat(l+1,:,1)=2/(T-l).*sum([m_It(1:end-l,:,:).*m_It_cycl(l+1:end,:,:)],1);
    conc_mat(l+1,:)=[(1/(T-l).*sum([It(1:end-l,:,1).*It(l+1:end,:,2)]+...
                 [(1-It(1:end-l,:,1)).*(1-It(l+1:end,:,2))],1)) ... 
                 (1/(G-l).*sum([It_g(1:end-l,:,1).*It_g(l+1:end,:,2)]+...
                 [(1-It_g(1:end-l,:,1)).*(1-It_g(l+1:end,:,2))],1))];
end

conc_matd=[reshape(conc_mat(1:sIt),sIt,1) reshape(conc_mat(sIt+1:end),2,sItg/2)'];
[strvcat(scen) num2str(conc_matd,'%5.2f %5.2f %5.2f')]





