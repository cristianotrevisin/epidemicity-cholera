%%% TO CREATE SCENARIOS WHERE MULTIPLE PORTIONS OF THE HAITIAN POPULATION IS PROGRESSIVELY VACCINATED
% TIME SPAN OF THE VACCINATION CAMPAIGN: 2 YEARS, STARTING ON NOV-01-2011 
% TIME BETWEEN 2 DOSES: 2 WEEKS (AS USUAL)
% 

close all
clear all
clc

[rSeq1, rSeq2] = load_data();
opt = "best";
[x1,x2] = get50(opt,rSeq1, rSeq2);

% run baseline model
[cases_AD1_week, time, y, cati, ocv] = SIARBV(2,2,x1,x2);
[Rt, et] = diagnosis(y, 2, 2, x1,x2);


et = smoothdata(et,'movmean',28);
Rt = smoothdata(Rt,'movmean',28);

% run model without vaccinations
[cases_AD1_week2, time2, y2, cati2, ocv2] = SIARBV(1,2, x1,x2);
[Rt2, et2] = diagnosis(y2,1,2,x1,x2);

et2 = smoothdata(et2,'movmean',28);
Rt2 = smoothdata(Rt2,'movmean',28);

% run model with 95% of the vaccinated population
[cases_AD1_week3, time3, y3, cati3, ocv3] = SIARBV(10,2,x1,x2);
[Rt3, et3] = diagnosis(y3,10,2,x1,x2);

et3 = smoothdata(et3,'movmean',28);
Rt3 = smoothdata(Rt3,'movmean',28);

% run model with 50% of the vaccinated population
[cases_AD1_week4, time4, y4, cati4, ocv4] = SIARBV(13,2,x1,x2);
[Rt4, et4] = diagnosis(y4,13,2,x1,x2);

et4 = smoothdata(et4,'movmean',28);
Rt4 = smoothdata(Rt4,'movmean',28);

% run model with 20% of the vaccinated population
[cases_AD1_week5, time5, y5, cati5, ocv5] = SIARBV(14,2,x1,x2);
[Rt5, et5] = diagnosis(y5,14,2,x1,x2);

et5 = smoothdata(et5,'movmean',28);
Rt5 = smoothdata(Rt5,'movmean',28);

%% FIND LAST NONZERO
timed = datenum('2010-10-20'):1:datenum('2017-07-01');

if length(et2)==3008
    timed2 = datenum('2010-10-20'):1:datenum('2019-01-13');
else 
    timed2 = datenum('2010-10-20'):1:datenum('2017-07-01');
end

cases_AD0_week = sum(cases_AD1_week,1);
l1 = find(cases_AD0_week>=1,1,'last'); m1 = find(timed==time(l1));
cases_AD0_week2 = sum(cases_AD1_week2,1);
l2 = find(cases_AD0_week2>=1,1,'last'); m2 = find(timed2==time(l2));
cases_AD0_week3 = sum(cases_AD1_week3,1);
l3 = find(cases_AD0_week3>=1,1,'last'); m3 = find(timed2==time(l3));
cases_AD0_week4 = sum(cases_AD1_week4,1);
l4 = find(cases_AD0_week4>=1,1,'last'); m4 = find(timed2==time(l4));
cases_AD0_week5 = sum(cases_AD1_week5,1);
l5 = find(cases_AD0_week5>=1,1,'last'); m5 = find(timed2==time(l5));
l0 = find(time==734791); m0 = find(timed==734791);

%% DEFINE COLOURS
c1 = 'k';
c2 = '#88b04b';
c3 = '#f4512c';
c4 = '#cfb095';
c5 = '#578ca9';

%% PLOT
cases_week = csvread('../data/cases.csv',1,1)'; 


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
        plot(time(1:l1),cases_AD0_week(1:l1),'-','color','black','linewidth',0.75); if l1 <350; plot(time(l1),cases_AD0_week(l1),'+k','LineWidth',2.5); end;
        plot(time2(52:l2),cases_AD0_week2(52:l2),'-','color',c2,'linewidth',0.75); if l2 <350; plot(time2(l2),cases_AD0_week2(l2),'+','color',c2,'LineWidth',2.5); end;
        plot(time3(52:l3),cases_AD0_week3(52:l3),'-','color',c3,'linewidth',0.75);if l3 <350; plot(time3(l3),cases_AD0_week3(l3),'+','color',c3,'LineWidth',2.5); end;
        plot(time4(52:l4),cases_AD0_week4(52:l4),'-','color',c4,'linewidth',0.75);if l4 <350; plot(time4(l4),cases_AD0_week4(l4),'+','color',c4,'LineWidth',2.5); end;
        plot(time5(52:l5),cases_AD0_week5(52:l5),'-','color',c5,'linewidth',0.75);if l5 <350; plot(time5(l5),cases_AD0_week5(l5),'+','color',c5,'LineWidth',2.5); end;
        ylabel('$C$', 'fontsize', 11, 'interpreter','latex')
        set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec,'Xticklabel',[])
        ylim([0 22500])
        box off
        hold off

    ax2 = axes('Position',[0.315 0.75 0.5 0.15]);
        text(0.95,0.9,'(b)','Units','normalized','FontSize',11)
        hold on
        b1 = stairs(timed,cumsum(sum(ocv.rv_1d,1) + sum(ocv.rv_2d,1)),'k','linewidth',0.75);
        b2 = stairs(timed2,cumsum(sum(ocv2.rv_1d,1) + sum(ocv2.rv_2d,1)),'color',c2,'linewidth',0.75);
        b3 = stairs(timed2,cumsum(sum(ocv3.rv_1d,1) + sum(ocv3.rv_2d,1)),'color',c3,'linewidth',0.75);
        b4 = stairs(timed2,cumsum(sum(ocv4.rv_1d,1) + sum(ocv4.rv_2d,1)),'color',c4,'linewidth',0.75);
        b5 = stairs(timed2,cumsum(sum(ocv5.rv_1d,1) + sum(ocv5.rv_2d,1)),'color',c5,'linewidth',0.75);
        ylabel('OCV (total)', 'fontsize',11)   
        l=legend([b1 b2 b3 b4 b5],'Baseline', 'No OCV', '95%', '50%', '20%','NumColumns',2);
        legend boxoff
        set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec)
        datetick('x','yyyy','keeplimits','keepticks')
        ylabel('OCVs (total)')
        hold off
        box off
    
    ax3 = axes('Position',[0.1 0.325 0.8 0.175]);
        text(0.95,0.9,'(c)','Units','normalized','FontSize',11)
        hold on
        plot(timed2,ones(length(timed2),1), 'color', 'red', 'linewidth',0.3)
        plot(timed(1:m1),Rt(1:m1),'-k','linewidth',0.75); if m1 ~= length(Rt); plot(timed(m1:end),Rt(m1:end),':k','linewidth',0.75); end
        plot(timed2(m0:m2),Rt2(m0:m2),'-','color',c2,'linewidth',0.75); if m2 ~= length(Rt2); plot(timed2(m2:end),Rt2(m2:end),':','color',c2,'linewidth',0.75); end
        plot(timed2(m0:m3),Rt3(m0:m3),'-','color',c3,'linewidth',0.75); if m3 ~= length(Rt3); plot(timed2(m3:end),Rt3(m3:end),':','color',c3,'linewidth',0.75); end
        plot(timed2(m0:m4),Rt4(m0:m4),'-','color',c4,'linewidth',0.75); if m4 ~= length(Rt4); plot(timed2(m4:end),Rt4(m4:end),':','color',c4,'linewidth',0.75); end
        plot(timed2(m0:m5),Rt5(m0:m5),'-','color',c5,'linewidth',0.75); if m5 ~= length(Rt5); plot(timed2(m5:end),Rt5(m5:end),':','color',c5,'linewidth',0.75); end
        ylabel('$\mathcal{R}_t$', 'fontsize',11, 'interpreter','latex')
        set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec,'Xticklabel',[])
        hold off
        box off
    
    ax4 = axes('Position',[0.1 0.1 0.8 0.175]);
        text(0.95,0.9,'(d)','Units','normalized','FontSize',11)
        hold on
        plot(timed2,zeros(length(timed2),1), 'color', 'red', 'linewidth',0.3)
        plot(timed(1:m1),et(1:m1),'-k','linewidth',0.75); if m1 ~= length(et); plot(timed(m1:end),et(m1:end),':k','linewidth',0.75); end
        plot(timed2(m0:m2),et2(m0:m2),'-','color',c2,'linewidth',0.75); if m2 ~= length(et2); plot(timed2(m2:end),et2(m2:end),':','color',c2,'linewidth',0.75); end
        plot(timed2(m0:m3),et3(m0:m3),'-','color',c3,'linewidth',0.75); if m3 ~= length(et3); plot(timed2(m3:end),et3(m3:end),':','color',c3,'linewidth',0.75); end
        plot(timed2(m0:m4),et4(m0:m4),'-','color',c4,'linewidth',0.75); if m4 ~= length(et4); plot(timed2(m4:end),et4(m4:end),':','color',c4,'linewidth',0.75); end
        plot(timed2(m0:m5),et5(m0:m5),'-','color',c5,'linewidth',0.75); if m5 ~= length(et5); plot(timed2(m5:end),et5(m5:end),':','color',c5,'linewidth',0.75); end
        ylabel('$e_t$ [1/d]', 'fontsize',11, 'interpreter','latex')  
        set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec)
        datetick('x','mmm-yy','keeplimits','keepticks')
        ylim([-0.2 0.8])
        hold off
        box off
        
    f.Units='points';
    f.Position=[0 0 400 400];