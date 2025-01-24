%program to forecast
fid=fopen('forecast_outsample.out','a');
randn('state',999);
rand('state',999);

freq=4; %observation frequency
Ni=1; %number of years to forecast
ser_fc={'YER';'KTAUF'}%;'PCR'yx;'ITR';'MTR';'XTR'}; %for which series do you want a forecast


fully_b=1; %flag to obtain fully Bayesian forecasts (taking into account estimation and parameter uncertainty)
stdep=1; %flag for state-conditional forecasts
H=0; %H to evaluate forecasts
Hf=freq*Ni; %forecast horizon 
Mf=1000; %number of forecasts

file='bds';
span=18;
Pool=1;
Shrink=0;
Sg_change=0;
Unit_spec_var=1;
perm=0;
ifig=7;
for cal_beg=[1988];
    cal_end=cal_beg+span;
    for cal_end_EST=[2000:0.25:2006.75];
    for KK=[2];
        for Lag_dlo=[2]; ['lag=' int2str(Lag_dlo)] %lag endogenous variable
            Lag_dir=0; %lag exogenous variable

            load('-mat',[file int2str(Pool) '_shr' int2str(Shrink) '_var' int2str(Unit_spec_var) '_perm' int2str(perm) '_sg' int2str(Sg_change) '_K' int2str(KK) 'end' int2str(Lag_dlo) '_ex' int2str(Lag_dir) '_grspecstd_' num2str(cal_end_EST,'%5.2f') '_' int2str(cal_end)]);
            %             load(['logit_test_' int2str(cal_beg) '_' int2str(cal_end) '_iden'])
            run_data_endswit_logit_forecast;
            %             run_model_test_endswit
            vars=find(sum(strcmp(ser_sel(ones(size(ser_fc,1),1),:),ser_fc(:,ones(1,size(ser_sel,2)))),1)); %index of forecasted series
            if H>0;
                cal_fore=[cal_fore;cal_fore(end)+[1:Hf]'./freq];
            elseif H==0
                cal_fore=cal_fore(end)+[1:Hf]'./freq;
            end
               
            [bfore_mat,bforeIMS_mat]=bayes_forecast(H,Hf,cal_fore,Mf,vars,fully_b,stdep,y,Z,alphamc,Qmc,sepsmc,lambmc,etaMSmc,IMSmc,ISmc,pool,unit_spec_var,K,dd,d,dMS,isp_lag_dlo);
            plot_forecast        
        end
    end
    end
end
fclose(fid)

