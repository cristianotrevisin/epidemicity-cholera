close all
clear all
ModelName = 'm4c';
load TEMP_RES.mat Sequences output

% close all
index0=find(squeeze(Sequences(:,end-1,1)),1,'first');
index1=find(squeeze(Sequences(:,1,1)),1,'last');

ribbon = 0;
Sequences =  Sequences(index0:index1,:,:);

[iter, chain] = find(squeeze(Sequences(:,end-1,:))==max(max(squeeze(Sequences(:,end-1,:)))));

vec = Sequences(max(iter),1:end-2,max(chain))

[cases_AD1_week, time2, y2] = m4c(vec);

time2 = datenum(2010,10,20):1:datenum(2010,10,20)+size(y2,1)-1;

figure()
subplot(4,1,1)
title('Susceptibles')
hold on
for i = 1:10
    plot(time2, y2(:,1 + (i-1)*5))
end
set(gca,'Xticklabel',[])
subplot(4,1,2)
title('Infected')
hold on
for i = 1:10
    plot(time2, y2(:,2 + (i-1)*5))
end
 set(gca,'Xticklabel',[])
subplot(4,1,3)
title('Recovered')
hold on
for i = 1:10
    plot(time2, y2(:,3 + (i-1)*5))
end
datetick('x','mmm-yy','keeplimits','keepticks')
subplot(4,1,4)
title('Bacteria')
hold on
for i = 1:10
    plot(time2, y2(:,4 + (i-1)*5))
end
datetick('x','mmm-yy','keeplimits','keepticks')
