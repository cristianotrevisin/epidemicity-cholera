%%% TO CREATE SCENARIOS WHERE THE NPIs ARE MULTIPLES OF THE ACTUAL NPIs
% 

close all
clear all
clc

[rSeq1, rSeq2] = load_data();
opt = "best";
[x1,x2] = get50(opt,rSeq1, rSeq2);

opt_ocv = 2;
% run baseline model
[cases_AD1_week, time, y, cati, ocv] = SIARBV(2,2,x1,x2);
[Rt, et] = diagnosis(y, 2, 2, x1,x2);


et = smoothdata(et,'movmean',28);
Rt = smoothdata(Rt,'movmean',28);

% run model without NPIs
[cases_AD1_week2, time2, y2, cati2, ocv2] = SIARBV(opt_ocv,1,x1,x2);
[Rt2, et2] = diagnosis(y2,opt_ocv,1, x1,x2);

et2 = smoothdata(et2,'movmean',28);
Rt2 = smoothdata(Rt2,'movmean',28);

% run model with twice as NPIs
[cases_AD1_week3, time3, y3, cati3, ocv3] = SIARBV(opt_ocv,3,x1,x2);
[Rt3, et3] = diagnosis(y3,opt_ocv,3, x1,x2);

et3 = smoothdata(et3,'movmean',28);
Rt3 = smoothdata(Rt3,'movmean',28);

% run model with 5 times as NPIs
[cases_AD1_week4, time4, y4, cati4, ocv4] = SIARBV(opt_ocv,9,x1,x2);
[Rt4, et4] = diagnosis(y4,opt_ocv,9, x1,x2);

et4 = smoothdata(et4,'movmean',28);
Rt4 = smoothdata(Rt4,'movmean',28);

%run model with 10 times as NPIs
[cases_AD1_week5, time5, y5, cati5, ocv5] = SIARBV(opt_ocv,4,x1,x2);
[Rt5, et5] = diagnosis(y5,opt_ocv,4, x1,x2);

et5 = smoothdata(et5,'movmean',28);
Rt5 = smoothdata(Rt5,'movmean',28);

%% PLOT
cases_week = csvread('../data/cases.csv',1,1)'; 
timed = datenum('2010-10-20'):1:datenum('2017-07-01');

if length(et2)==3008
    timed2 = datenum('2010-10-20'):1:datenum('2019-01-13');
else 
    timed2 = datenum('2010-10-20'):1:datenum('2017-07-01');
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
        plot(time,sum(cases_AD1_week),'-k','linewidth',0.75);
        plot(time2,sum(cases_AD1_week2),'-b','linewidth',0.75);
        plot(time3,sum(cases_AD1_week3),'-r','linewidth',0.75);
        plot(time4,sum(cases_AD1_week4),'-g','linewidth',0.75);
        plot(time5,sum(cases_AD1_week5),'-m','linewidth',0.75);
        ylabel('$\Delta C$', 'fontsize', 11, 'interpreter','latex')
        set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec,'Xticklabel',[])
        ylim([0 22500])
        box off
        hold off

    ax2 = axes('Position',[0.315 0.75 0.5 0.15]);
        text(0.95,0.9,'(b)','Units','normalized','FontSize',11)
        hold on
        b5 = bar(timed,sum(diff([zeros(1,10) ; cati5]),2),'m');
        b5.FaceAlpha = 0.75;
        b4 = bar(timed,sum(diff([zeros(1,10) ; cati4]),2),'g');
        b4.FaceAlpha = 0.75;
        b3 =  bar(timed,sum(diff([zeros(1,10) ; cati3]),2),'r');
        b3.FaceAlpha = 0.75;
        b2 = bar(timed2,sum(diff([zeros(1,10) ; cati2]),2),'b');
        b2.FaceAlpha = 0.75;
        b1 = bar(timed,sum(diff([zeros(1,10) ; cati]),2),'k');
        b1.FaceAlpha = 0.75;
        l=legend([b1 b2 b3 b4 b5], {'Baseline', 'No NPI', 'NPI x 2', 'NPI x 5', 'NPI x 10'}, 'location', 'northwest');
        h=findobj(l,'type','patch'); 
        set(h,'ydata',[0.,0.5,0.5,0.3,0.3],'xdata',[0.2,0.2,0.4,0.4,0.2]);
        legend boxoff
        set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec)    
        datetick('x','yyyy','keeplimits','keepticks')
        ylabel('Weekly NPI')
        hold off
        box off
    
    ax3 = axes('Position',[0.1 0.325 0.8 0.175]);
        text(0.95,0.9,'(c)','Units','normalized','FontSize',11)
        hold on
        plot(timed2,ones(length(timed2),1), 'color', 'red', 'linewidth',0.3)
        plot(timed,Rt,'-k','linewidth',0.75)
        plot(timed2,Rt2,'-b','linewidth',0.75)
        plot(timed,Rt3,'-r','linewidth',0.75)
        plot(timed,Rt4,'-g','linewidth',0.75)
        plot(timed,Rt5,'-m','linewidth',0.75)
        ylabel('$\mathcal{R}_t$', 'fontsize',11, 'interpreter','latex')
        set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec,'Xticklabel',[])
        hold off
        box off
    
    ax4 = axes('Position',[0.1 0.1 0.8 0.175]);
        text(0.95,0.9,'(d)','Units','normalized','FontSize',11)
        hold on
        plot(timed2,zeros(length(timed2),1), 'color', 'red', 'linewidth',0.3)
        plot(timed,et,'-k','linewidth',0.75)
        plot(timed2,et2,'b','linewidth',0.75)
        plot(timed,et3,'r','linewidth',0.75)
        plot(timed,et4,'g','linewidth',0.75)
        plot(timed,et5,'m','linewidth',0.75)
        ylabel('$e_t$ [1/d]', 'fontsize',11, 'interpreter','latex')  
        set(gca,'Xlim',[timed(1) timed2(end)],'Xtick',tick_vec)
        datetick('x','mmm-yy','keeplimits','keepticks')
        ylim([-0.2 0.8])
        hold off
        box off

    f.Units='points';
    f.Position=[0 0 400 400];