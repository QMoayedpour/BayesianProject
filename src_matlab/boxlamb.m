ctries=ser_sel;
index=[1:1:size(lambmc,2)];
indexm=[1:5:size(lambmc,1)];
figure(gcf+1)
boxplot(lambmc(indexm,index),0,'+',[],5)
ylabel('\lambda_i','FontSize',12);
xlabel(' ');
ax1=gca;
set(ax1,'XLim',[0 length(index)+1],'XTickLabel',ctries)

figure(gcf+1)
boxplot([sepsmc(indexm,ones(1,length(index)))./lambmc(indexm,index)],0,'+',[],5)
ylabel('\sigma_i^2','FontSize',12);
xlabel(' ');
ax1=gca;
set(ax1,'XLim',[0 length(index)+1],'XTickLabel',ctries)
