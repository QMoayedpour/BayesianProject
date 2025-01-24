conf_ind=((cat_dum==2)+(cat_dum==3))==1;
ser_lead=ser_sel(conf_ind);
lead_ind=find(conf_ind);
pl_ind=lead_ind(find(prob(2,conf_ind)));
pc_ind=prob(1,:)==1;

rec=2; %1 means below-average growth, 2 means recession
nst=2;
it=size(LIMSmc,2);
KIMS=size(LIMSmc,1);

%probMS=zeros(nst,it,K);
probMS=1-mean(LIMSmc,3);
if unit_spec_var
    y_gr=dlo(end-it+1:end,:);%./ kron(ones(it,1),sqrt(mean(sepsmc(:,ones(1,size(lambmc,2)))./lambmc,1)));
    data_raw=data_raw(end-it+1:end,:);
    data_raw(:,conf_ind)=(data_raw(:,conf_ind)-kron(ones(it,1),mean(data_raw(:,conf_ind),1)))./ kron(ones(it,1),std(data_raw(:,conf_ind),[],1));
elseif ~unit_spec_var
    y_gr=dlo(end-it+1:end,:)./ kron(ones(it,1),sqrt(mean(sepsmc(:,ones(1,size(lambmc,2))),1)));
end

time = 1:it;
c=cal_end_est+1/4;
cal = sort(c-time/4)';


figure(gcf+1)
for k=2
%     subplot(3,1,1)
% 
%     bar(cal,probMS(1,:),0.2);
%     hold on
%     hold on;
%     if ~DStruc
%         ylabel({'Coincident series';'I_1_t=1'},'FontSize',12)
%         %ylabel(['Differences' ],'FontSize',12)
%     elseif DStruc
%         ylabel(['\rho=' num2str(mean(rhomc(:,k)),'%f5.2')])
%     end
%     ax1=gca;
%     set(ax1,'XLim',[cal(1) cal(end)],'XTick',[cal_beg:2:cal_end],'YLim',[0 1],'YTick',[0:0.5:1])
%     ax2=axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none');
%     set(ax2,'XLim',[cal(1) cal(end)],'XTickLabel',' ','YLim',...
%         [min(min(y_gr(:,pc_ind)))-0.5 max(max(y_gr(:,pc_ind)))+0.5])
% 
%     if any(prob(k,:))
%         line(cal,y_gr(:,pc_ind),'Parent',ax2,'Color','k')
%         hold on
%         line(cal,y_gr(:,1),'Parent',ax2,'Color','r')
%         %
%         %         hold on
%     end
 

    subplot(2,1,1)
    bar(cal,probMS(1,:),0.2);
    hold on
    hold on;
    if ~DStruc
        ylabel({'Differenced';'I_1_t=1'},'FontSize',12)
        %ylabel(['Differences' ],'FontSize',12)
    elseif DStruc
        ylabel(['\rho=' num2str(mean(rhomc(:,k)),'%f5.2')])
    end
    ax1=gca;
    set(ax1,'XLim',[cal(1) cal(end)],'XTick',[cal_beg:2:cal_end],'YLim',[0 1],'YTick',[0:0.5:1])
    ax2=axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none');
    set(ax2,'XLim',[cal(1) cal(end)],'XTickLabel',' ','YLim',...
        [min(min(y_gr(:,pl_ind)))-0.5 max(max(y_gr(:,pl_ind)))+0.5])

    if any(prob(k,:))
        line(cal,y_gr(:,pl_ind),'Parent',ax2,'Color','k')
%         hold on
%         line(cal,y_gr(:,1),'Parent',ax2,'Color','r')
        %
        %         hold on
    end
    %     line(cal,zeros(length(cal),1),'Parent',ax2,'Color','k')
    %    line(cal,ones(length(cal),1)*mean(mean(y_gr(:,prob(k,:)==1))),'Parent',ax2,'Color','k')

    %     if k==2 title(['I_t=1']);end
    subplot(2,1,2)

    bar(cal,probMS(2,:),0.2);
    hold on
    if ~DStruc
        ylabel({'Levels';'I_2_t=1'},'FontSize',12)
    elseif DStruc
        ylabel(['\rho=' num2str(mean(rhomc(:,k)),'%f5.2')])
    end
    ax1=gca;
    set(ax1,'XLim',[cal(1) cal(end)],'XTick',[cal_beg:2:cal_end],'YLim',[0 1],'YTick',[0:0.5:1])
    ax2=axes('Position',get(ax1,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none');
    set(ax2,'XLim',[cal(1) cal(end)],'XTickLabel',' ','YLim',...
        [min(min(data_raw(:,pl_ind)))-0.5 max(max(data_raw(:,pl_ind)))+0.5])
    %     %     set(ax2,'XLim',[cal(1) cal(end)],'XTickLabel',' ','YLim',[-5 5])
    %     %     if k==1 set(ax2,'XLim',[cal(1) cal(end)],'YLim',[-5 5]); end
    %     %     if find(y_gr(:,prob(k,:)==1))
    if any(prob(k,:))
        line(cal,data_raw(:,pl_ind),'Parent',ax2,'Color','k')
        %         hold on
        %         line(cal,y_gr(:,1),'Parent',ax2,'Color','r')
        %
        %         hold on
    end


end


prob_lead=prob(:,cat_dum==2);
