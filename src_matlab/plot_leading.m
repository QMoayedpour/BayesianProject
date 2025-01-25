%program to plot some interesting output for various estimated models

% cal_beg=[1995];
% cal_end=[1986:2003];
span=18;
pool=1;
poolparms=1;
shrink=1;
sg_change=1;
unit_spec_var=1;
perm=0;
File='logit';
DSTruc=0;
S_Logit=1;
Pclas=1;
rest_diff=0;
for cal_beg=[1988];
    cal_end=cal_beg+span;
    for cal_end_Est=[2006.75];
    for K=[3];
        for lag_dlo=[2]; ['lag=' int2str(lag_dlo)] %lag endogenous variable
            lag_dir=0; %lag exogenous variable
            
            %             load('-mat',['logit' int2str(pool) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_est,'%5.2f') '_' int2str(cal_end) '_iden']);
%             if Pclas
%                 load('-mat',[File int2str(S_Logit) '_dyn' int2str(DSTruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_Est,'%5.2f') '_' int2str(cal_end) '_iden']);
%             elseif ~Pclas
                load('-mat',[File int2str(S_Logit) '_dyn' int2str(DSTruc) '_shr' int2str(shrink) '_var' int2str(unit_spec_var) '_perm' int2str(perm) '_sg' int2str(sg_change) '_K' int2str(K) 'end' int2str(lag_dlo) '_ex' int2str(lag_dir) '_grspecstd_' num2str(cal_end_Est,'%5.2f') '_' int2str(cal_end) '_iden']);
%             end
                %             load(['logit_test_' int2str(cal_beg) '_' int2str(cal_end) '_iden'])
            eval=1;
            Q=Q0;
            run_model4_endswit_logit;
%             run_model_test_endswit
            plotprob;
            
            plot_leadser
           clear functions
            % ctr=ctr+1;
        end
    end
    end
end

