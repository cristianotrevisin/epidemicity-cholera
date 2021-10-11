close all
clear all
clc

% 1 SCENARIO 0 -> NO VACCINATION
% 2 SCENARIO B -> WHAT ACTUALLY HAPPENED
% 3 SCENARIO 1 -> TWO DEPARTMENTS (CT, AT) 3.71 M DOSES OVER 730 DAYS AU PRORATA
% 4 SCENARIO 2 -> THREE DEPTS (CT, AT, OU) 9.76 M DOSES OVER 730 DAYS AU PRORATA
% 5 SCENARIO 3 -> SLOW NATIONAL 16.37 M DOSES OVER 1825 DAYS AU PRORATA
% 6 SCENARIO 4 -> FAST NATIONAL 16.37 M DOSES OVER 730 DAYS AU PRORATA
% 7 SCENARIO 5 -> HIGH COVERAGE 20.91 M DOSES OVER 730 DAYS AU PRORATA
% 8 SCENARIO b1-> TWICE AS MANY VACCINATIONS
% 9 SCENARIO b2-> VACCINATIONS STARTING ON NOV 1st

% SECOND DOSE ALWAYS GIVEN TWO WEEK AFTER THE FIRST DOSE

% SCENARIOS 1-4: 70% COVERAGE WITH TWO DOSES; 10% COVERAGE WITH ONE DOSE ->
% OF ALL DOSES, 80/150 ARE FIRST DOSES; 70/150 ARE SECOND DOSES

% SCENARIO 5: 95% COVERAGE WITH TWO DOSES; 1.67% COVERAGE WITH ONE DOSE ->
% OF ALL DOSES, 96.67/191.67 ARE FIRST DOSES; 95/191.67 ARE SECOND DOSES


%1 SCENARIO 0 -> NO NPI
%2 SCENARIO B -> WHAT ACTUALLY HAPPENED
%3 SCENARIO 1 -> TWICE AS NPI
%4 SCENARIO 2 -> TENFOLD AS NPI
%5 SCENARIO 3 -> 1 NPI EVERY 100 CASES
%6 SCENARIO 4 -> 1 NPI EVERY 50 CASES
%7 SCENARIO 5 -> 1 NPI EVERY 10 CASES
%8 SCENARIO 1b -> ANTICIPATED NPI

opt_ocv = 2;
opt_npi = 7;

[cases_AD1_week, time, y, cati, ocv] = m4s(2,2);
[Rt, et] = diagnosis_MOD(y, 2, 2);


et = smoothdata(et,'movmean',28);
Rt = smoothdata(Rt,'movmean',28);

[cases_AD1_week2, time2, y2, cati2, ocv2] = m4s(opt_ocv,opt_npi);
[Rt2, et2] = diagnosis_MOD1(y2,opt_ocv,opt_npi,cati2);

et2 = smoothdata(et2,'movmean',28);
Rt2 = smoothdata(Rt2,'movmean',28);

%% PLOT
close all
cases_week = csvread('data/cases.csv',1,1)'; 

timed = datenum('2010-10-20'):1:datenum('2017-07-01');

timed2 = datenum('2010-10-20'):1:datenum('2017-07-01');

if length(et2)==3008
    timed3 = datenum('2010-10-20'):1:datenum('2019-01-13');
else 
    timed3 = datenum('2010-10-20'):1:datenum('2017-07-01');
end
tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.11.2011','dd.mm.yyyy') ...
    datenum('01.11.2012','dd.mm.yyyy') datenum('01.11.2013','dd.mm.yyyy') ...
    datenum('01.11.2014','dd.mm.yyyy') datenum('01.11.2015','dd.mm.yyyy')...
    datenum('01.11.2016','dd.mm.yyyy') datenum('01.11.2017','dd.mm.yyyy')...
    datenum('01.11.2018','dd.mm.yyyy') datenum('01.11.2019','dd.mm.yyyy')];

tickvec = [-0.05, 0, 1, 100];
f = figure(4002);

    subplot(5,1,1) 
    text(0.95,0.9,'(a)','Units','normalized','FontSize',11)
    hold on
    set(gcf,'color','white')
    %p2 = bar(time, sum(cases_week(:,length(time))),'FaceColor','#C1C1C1')
    plot(time,sum(cases_AD1_week),'-k');
    plot(time2,sum(cases_AD1_week2),'--b');
    box off
    legend('Scenario 0', 'Scenario 1', 'location', 'north')
    legend boxoff
    %p2 = plot(time,sum(cases_week),'o','markeredgecolor','red','markerfacecolor','white','markersize',2,'color',[0.5 0.5 0.5],'linewidth', 1);
    %legend('Observation', 'Simulation', 'Scenario', 'fontsize', 12, 'location','east')
    ylabel('$\Delta C$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[time2(1) time2(end)],'Xtick',tick_vec)
    set(gca,'Xticklabel',[])
   % legend boxoff
    box off
    ylim([0 26000])
    hold off
%legend('Scenario 0', 'Doubled OCV', 'Anticipated OCV', 'Doubled WaSH', 'Anticipated WaSH', 'location', 'bestoutside')

    subplot(5,1,2)
    text(0.95,0.9,'(b)','Units','normalized','FontSize',11)
    hold on
    set(gca,'Xlim',[timed(1) timed(end)])   
    plot(timed,ones(length(timed),1), 'color', 'red', 'linewidth',0.3)
    plot(timed,Rt,'-k')
    plot(timed3,Rt2,'--b')
    %ylim([0.01 100])
    ylabel('$\mathcal{R}_t$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec) 
    set(gca,'Xticklabel',[])
    hold off
    box off
    
    subplot(5,1,3)
    text(0.95,0.9,'(c)','Units','normalized','FontSize',11)
    hold on
    set(gca,'Xlim',[timed(1) timed(end)])
        plot(timed,et,'-k')
        plot(timed,zeros(length(timed),1), 'color', 'red', 'linewidth',0.3)
        plot(timed3,et2,'--b')
        ylabel('$e_t$', 'fontsize',11, 'interpreter','latex')     
    ylim([-0.2 0.8])
set(gca,'Xticklabel',[])
box off

    subplot(5,1,4)
    text(0.95,0.9,'(d)','Units','normalized','FontSize',11)
    hold on
    set(gca,'Xlim',[timed(1) timed(end)])
    b1 = bar(timed,sum(diff([zeros(1,10) ; cati]),2),'k')
       b1.FaceAlpha = 0.75
       b2 = bar(timed3,sum(diff([zeros(1,10) ; cati2]),2),'b')
       b2.FaceAlpha = 0.75
        ylabel('NPI', 'fontsize',11)     
set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec) 
set(gca,'Xticklabel',[])
datetick('x','mmm-yy','keeplimits','keepticks')
box off

    subplot(5,1,5)
    text(0.95,0.9,'(e)','Units','normalized','FontSize',11)
    hold on
    set(gca,'Xlim',[timed(1) timed(end)])
%     bar(timed,sum(ocv.rv_1d,1) + sum(ocv.rv_2d,1),'k')
    stairs(timed,cumsum(sum(ocv.rv_1d,1) + sum(ocv.rv_2d,1)),'k')
        stairs(timed3,cumsum(sum(ocv2.rv_1d,1) + sum(ocv2.rv_2d,1)),'b')
        ylabel('OCV (doses)', 'fontsize',11)     
set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec) 
datetick('x','mmm-yy','keeplimits','keepticks')
box off
%line([timed(1) timed(end)], [0 0], 'color', 'magenta','LineWidth',0.1,'LineStyle','--')
hold off
f.Units='points';
f.Position=[0 0 400 600];





%%%%%%%%HOUSEKEEPING
% subplot(4,1,3)
%     text(0.95,0.9,'(c)','Units','normalized','FontSize',11)
%     hold on
%     Rtp = Rt; Rtp(Rtp<1) = NaN; Rtn = Rt; Rtn(Rtn>1) = NaN;
%     Rtp2 = Rt2; Rtp2(Rtp2<1) = NaN; Rtn2 = Rt2; Rtn2(Rtn2>1) = NaN;
%     set(gca,'Xlim',[timed(1) timed(end)])   
%     fill([timed(1) timed2(end) timed2(end) timed(1)], [0.01 0.01 1 1],'b','facealpha',0.075,'EdgeColor','none')
%     hold on
%     fill([timed(1) timed2(end) timed2(end) timed(1)], [1 1 100 100],'r','facealpha',0.075,'EdgeColor','none')
%     plot(timed,Rtp,'-b','Linewidth',1)
%     plot(timed,Rtn,'--b','Linewidth',1)
%     plot(timed2,Rtp2,'-r','Linewidth',1)
%     plot(timed2,Rtn2,'--r','Linewidth',1)
%     %ylim([0.01 100])
%     ylim([min(min(0.5./Rt), min(0.5./Rt2)) max(max(2*Rt), max(2*Rt2))])
%     set(gca, 'yscale', 'log')
%     yticks([0.2 1 5])
%     ylabel('$R_t$', 'fontsize',11, 'interpreter','latex')
%     set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec) 
%     set(gca,'Xticklabel',[])
%     hold off
%     box off
%     
%     subplot(4,1,4)
%     etp = et; etp(etp<0) = NaN; etn = et; etn(etn>0) = NaN;
%     etp2 = et2; etp2(etp2<0) = NaN; etn2 = et2; etn2(etn2>0) = NaN;
%     text(0.95,0.9,'(d)','Units','normalized','FontSize',11)
%     hold on
%     set(gca,'Xlim',[timed(1) timed(end)])
%     fill([timed2(1) timed2(end) timed2(end) timed2(1)], [-1 -1 0 0],'b','facealpha',0.075,'EdgeColor','none')
%     hold on
%     fill([timed2(1) timed2(end) timed2(end) timed2(1)], [0 0 max(max(log(1+et))+1, max(log(1+et2)+1)) max(max(log(1+et))+1, max(log(1+et2)+1))],'r','facealpha',0.075,'EdgeColor','none')
%    % ylim([-1 max(log(1+et))+1])
%         plot(timed,log(1+etp),'-b','Linewidth',1)
%         plot(timed,log(1+etn),'--b','Linewidth',1)
%         plot(timed2,log(1+etp2),'-r','Linewidth',1)
%         plot(timed2,log(1+etn2),'--r','Linewidth',1)
%         ylabel('$\log(1+e_t)$', 'fontsize',11, 'interpreter','latex')     
%     ylim([-1.2*max(max(abs(log(1+et))),max(abs(log(1+et2))))  1.2*max(max(abs(log(1+et))),max(abs(log(1+et2))))])
% set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec)    
% datetick('x','mmm-yy','keeplimits','keepticks')
% box off
% %line([timed(1) timed(end)], [0 0], 'color', 'magenta','LineWidth',0.1,'LineStyle','--')
% hold off
% f.Units='points';
% f.Position=[0 0 450 600];