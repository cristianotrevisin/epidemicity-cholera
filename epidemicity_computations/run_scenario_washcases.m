close all

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
%9 SCENARIO 1.5 -> 5FOLD AS NPI

opt_ocv = 2;
opt_npi = 7;

[cases_AD1_week, time, y, cati, ocv] = m4s(2,2);
[Rt, et] = diagnosis_MOD(y, 2, 2);


et = smoothdata(et,'movmean',28);
Rt = smoothdata(Rt,'movmean',28);

[cases_AD1_week2, time2, y2, cati2, ocv2] = m4s(opt_ocv,1);
[Rt2, et2] = diagnosis_MOD1(y2,opt_ocv,1, cati2);

et2 = smoothdata(et2,'movmean',28);
Rt2 = smoothdata(Rt2,'movmean',28);

[cases_AD1_week3, time3, y3, cati3, ocv3] = m4s(opt_ocv,5);
[Rt3, et3] = diagnosis_MOD1(y3,opt_ocv,5, cati3);

et3 = smoothdata(et3,'movmean',28);
Rt3 = smoothdata(Rt3,'movmean',28);

[cases_AD1_week4, time4, y4, cati4, ocv4] = m4s(opt_ocv,6);
[Rt4, et4] = diagnosis_MOD1(y4,opt_ocv,6, cati4);

et4 = smoothdata(et4,'movmean',28);
Rt4 = smoothdata(Rt4,'movmean',28);

[cases_AD1_week5, time5, y5, cati5, ocv5] = m4s(opt_ocv,7);
[Rt5, et5] = diagnosis_MOD1(y5,opt_ocv,7, cati5);

et5 = smoothdata(et5,'movmean',28);
Rt5 = smoothdata(Rt5,'movmean',28);

%%
close all
cases_week = csvread('data/cases.csv',1,1)'; 

timed = datenum('2010-10-20'):1:datenum('2017-07-01');

timed2 = datenum('2010-10-20'):1:datenum('2017-07-01');

a = datenum(2013,03,01);
b = datenum(2019,01,13);

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

    ax1 = axes('Position',[0.1 0.55 0.8 0.35]);
    text(0.95,0.9,'(a)','Units','normalized','FontSize',11)
    hold on
    set(gcf,'color','white')
    %p2 = bar(time, sum(cases_week(:,length(time))),'FaceColor','#C1C1C1')
    plot(time,sum(cases_AD1_week),'-k','linewidth',0.75);
    plot(time2,sum(cases_AD1_week2),'-b','linewidth',0.75);
    plot(time3,sum(cases_AD1_week3),'-r','linewidth',0.75);
    plot(time4,sum(cases_AD1_week4),'-g','linewidth',0.75);
    plot(time5,sum(cases_AD1_week5),'-m','linewidth',0.75);
    box off
    
    %p2 = plot(time,sum(cases_week),'o','markeredgecolor','red','markerfacecolor','white','markersize',2,'color',[0.5 0.5 0.5],'linewidth', 1);
    %legend('Observation', 'Simulation', 'Scenario', 'fontsize', 12, 'location','east')
    ylabel('$\Delta C$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[a b],'Xtick',tick_vec)
    set(gca,'Xticklabel',[])
   % legend boxoff
    box off
    ylim([0 22500])
    hold off
%legend('Scenario 0', 'Doubled OCV', 'Anticipated OCV', 'Doubled WaSH', 'Anticipated WaSH', 'location', 'bestoutside')

    ax2 = axes('Position',[0.315 0.75 0.5 0.15]);
    text(0.95,0.9,'(b)','Units','normalized','FontSize',11)
    hold on
    set(gca,'Xlim',[timed(1) timed3(end)])
    b5 = bar(timed3,sum(diff([zeros(1,10) ; cati5]),2),'m')
       b5.FaceAlpha = 0.75
       b4 = bar(timed3,sum(diff([zeros(1,10) ; cati4]),2),'g')
       b4.FaceAlpha = 0.75
       b3 =  bar(timed3,sum(diff([zeros(1,10) ; cati3]),2),'r')
       b3.FaceAlpha = 0.75
       b2 = bar(timed3,sum(diff([zeros(1,10) ; cati2]),2),'b')
       b2.FaceAlpha = 0.75
       b1 = bar(timed,sum(diff([zeros(1,10) ; cati]),2),'k')
       b1.FaceAlpha = 0.75
        ylabel('NPI', 'fontsize',11) 
       l=legend([b1 b2 b3 b4 b5],'  Baseline', '  No NPI', '  x 0.05', '  x 0.10', '  x 0.20', 'location', 'northwest')
       h=findobj(l,'type','patch'); 
set(h,'ydata',[0.,0.5,0.5,0.3,0.3],'xdata',[0.2,0.2,0.4,0.4,0.2]);
legend boxoff
set(gca,'Xlim',[timed(1)-100 timed3(end)],'Xtick',tick_vec) 
set(gca,'Xticklabel',[])
    ylabel('Weekly NPI')
    set(gca,'Xlim',[timed(1)-199 timed(end)])   
    datetick('x','yyyy','keeplimits','keepticks')
    hold off
    box off
    
    ax3 = axes('Position',[0.1 0.325 0.8 0.175]);
    text(0.95,0.9,'(c)','Units','normalized','FontSize',11)
    hold on
    set(gca,'Xlim',[timed(1) timed3(end)])   
    plot(timed3,ones(length(timed3),1), 'color', 'red', 'linewidth',0.3)
    plot(timed,Rt,'-k','linewidth',0.75)
    plot(timed3,Rt2,'-b','linewidth',0.75)
    plot(timed3,Rt3,'-r','linewidth',0.75)
    plot(timed3,Rt4,'-g','linewidth',0.75)
    plot(timed3,Rt5,'-m','linewidth',0.75)
    %ylim([0.01 100])
    ylabel('$\mathcal{R}_t$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[a b],'Xtick',tick_vec) 
    set(gca,'Xticklabel',[])
%     legend('Baseline', 'No NPI', 'NPI x 2', 'NPI x 5', 'NPI x 10', 'location', 'north')
%     legend boxoff
    hold off
    box off
    
    ax4 = axes('Position',[0.1 0.1 0.8 0.175]);
    text(0.95,0.9,'(d)','Units','normalized','FontSize',11)
    hold on
    set(gca,'Xlim',[timed(1) timed3(end)])
        plot(timed,et,'-k','linewidth',0.75)
        plot(timed3,zeros(length(timed3),1), 'color', 'red', 'linewidth',0.3)
        plot(timed3,et2,'b','linewidth',0.75)
        plot(timed3,et3,'r','linewidth',0.75)
        plot(timed3,et4,'g','linewidth',0.75)
        plot(timed3,et5,'m','linewidth',0.75)
        ylabel('$e_t$ [1/d]', 'fontsize',11, 'interpreter','latex')  
        set(gca,'Xlim',[a b],'Xtick',tick_vec)
        datetick('x','mmm-yy','keeplimits','keepticks')
    ylim([-0.2 0.8])

box off

%     subplot(5,1,4)
%     text(0.95,0.9,'(d)','Units','normalized','FontSize',11)
%     hold on
%     set(gca,'Xlim',[timed(1) timed3(end)])
%     b1 = bar(timed,sum(diff([zeros(1,10) ; cati]),2),'k')
%        b1.FaceAlpha = 0.75
%        b2 = bar(timed3,sum(diff([zeros(1,10) ; cati2]),2),'b')
%        b2.FaceAlpha = 0.75
%        b3 =  bar(timed2,sum(diff([zeros(1,10) ; cati3]),2),'r')
%        b3.FaceAlpha = 0.75
%        b4 = bar(timed2,sum(diff([zeros(1,10) ; cati4]),2),'g')
%        b4.FaceAlpha = 0.75
%        b8 = bar(timed2,sum(diff([zeros(1,10) ; cati5]),2),'m')
%        b8.FaceAlpha = 0.75
%         ylabel('NPI', 'fontsize',11)     
% set(gca,'Xlim',[a b],'Xtick',tick_vec) 
% set(gca,'Xticklabel',[])
% datetick('x','mmm-yy','keeplimits','keepticks')
% box off
% 
%     subplot(5,1,5)
%     text(0.95,0.9,'(e)','Units','normalized','FontSize',11)
%     hold on
%     set(gca,'Xlim',[timed(1) timed3(end)])
% %     bar(timed,sum(ocv.rv_1d,1) + sum(ocv.rv_2d,1),'k')
%     stairs(timed,cumsum(sum(ocv.rv_1d,1) + sum(ocv.rv_2d,1)),'k')
%         stairs(timed3,cumsum(sum(ocv2.rv_1d,1) + sum(ocv2.rv_2d,1)),'b')
%         stairs(timed2,cumsum(sum(ocv3.rv_1d,1) + sum(ocv3.rv_2d,1)),'r')
%         stairs(timed2,cumsum(sum(ocv4.rv_1d,1) + sum(ocv4.rv_2d,1)),'g')
%         stairs(timed2,cumsum(sum(ocv5.rv_1d,1) + sum(ocv5.rv_2d,1)),'m')
%         ylabel('OCV (doses)', 'fontsize',11)     
% set(gca,'Xlim',[a b],'Xtick',tick_vec) 
% datetick('x','mmm-yy','keeplimits','keepticks')
% box off
% %line([timed(1) timed(end)], [0 0], 'color', 'magenta','LineWidth',0.1,'LineStyle','--')
% hold off
f.Units='points';
f.Position=[0 0 400 400];





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