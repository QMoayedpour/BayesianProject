%program to forecast
fid=fopen('forecast.out','w');
randn('state',999);
rand('state',999);

freq=4; %observation frequency
Ni=1; %number of years to forecast
ser_fc={'YER';'KTAUF'}%;'PCR';'ITR';'MTR';'XTR'}; %for which series do you want a forecast
end_for=1; %flag to compute out-of-sample probabilistic forecasts at the end of the observation period


fully_b=0; %flag to obtain fully Bayesian forecasts (taking into account estimation and parameter uncertainty)
stdep=1; %flag for state-conditional forecasts
H=24; %H to evaluate forecasts
Hf=freq*Ni; %forecast horizon 
Mf=1000; %number of forecasts

span=18;
pool=1;
Shrink=1;
Unit_spec_var=1;
sg_change=1;
File='logit';
DSTruc=0;
S_Logit=1;
rest_diff=0;
perm=0;
ifig=7;
for cal_beg=[1988];
    cal_end=cal_beg+span;
    for cal_end_est=[2006.75];
    for K=[3];
        for lag_dlo=[2]; ['lag=' int2str(lag_dlo)] %lag endogenous variable
            lag_dir=0; %lag exogenous variable

            load('-mat',['logit' int2str(pool) '_dyn' int2str(DSTruc) '_shr' int2str(Shrink) '_var' int2str(Unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_est,'%5.2f') '_' int2str(cal_end)]);
            %             load(['logit_test_' int2str(cal_beg) '_' int2str(cal_end) '_iden'])
            run_data_endswit_logit_forecast;
            %             run_model_test_endswit
            vars=find(sum(strcmp(ser_sel(ones(size(ser_fc,1),1),:),ser_fc(:,ones(1,size(ser_sel,2)))),1)); %index of forecasted series
            cal_fore=[cal_fore;cal_fore(end)+[1:Hf]'./freq];
            [bfore_mat,bforeIMS_mat]=bayes_forecast(H,Hf,cal_fore,Mf,vars,fully_b,stdep,y,Z,alphamc,Qmc,sepsmc,lambmc,etaMSmc,IMSmc,ISmc,pool,unit_spec_var,K,dd,d,dMS,isp_lag_dlo);
            plot_forecast
            if end_for
                prob_forecasts
            end
        end
    end
    end
end
fclose(fid);

