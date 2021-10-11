close all
clear all
clc

load ../data/geodata POPnodes WS_dept
POPnodes = POPnodes' * WS_dept;
POPnodes = POPnodes';



cases_week = csvread('../data/cases.csv',1,1)';
cases_week = cases_week(:,1:350); 
wm=10;
[cases_AD1_week, time, y] = m5c_washmult(1);
[Rt, et] = diagnosis(y,1,1);



Rt = smoothdata(Rt,'movmean',28);
et = smoothdata(et,'movmean',28);

%% PLOT
close all
timed = datenum('2010-10-20'):1:datenum('2017-07-01');
tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.11.2011','dd.mm.yyyy') ...
    datenum('01.11.2012','dd.mm.yyyy') datenum('01.11.2013','dd.mm.yyyy') ...
    datenum('01.11.2014','dd.mm.yyyy') datenum('01.11.2015','dd.mm.yyyy')...
    datenum('01.11.2016','dd.mm.yyyy')];

tickvec = [-0.05, 0, 1, 100];

indexswitchRt = [1 find(diff(sign(Rt-1)))+1];
indexswitchet = [1 find(diff(sign(et)))+1];


f = figure(4002);

    subplot(3,1,1) 
    text(0.95,0.9,'(a)','Units','normalized','FontSize',11)
    hold on
    set(gcf,'color','white')
    p2 = bar(time, sum(cases_week),'FaceColor','#C1C1C1')
    plot(time,sum(cases_AD1_week),'color','black','Linewidth',1);
    box off
    %p2 = plot(time,sum(cases_week),'o','markeredgecolor','red','markerfacecolor','white','markersize',2,'color',[0.5 0.5 0.5],'linewidth', 1);
    %legend('Observation', 'Simulation','fontsize', 11, 'location','southeast')
    ylabel('$\Delta C$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    set(gca,'Xticklabel',[])
    %legend boxoff
    box off
    ylim([0 26000])
    hold off

    subplot(3,1,2)
    text(0.95,0.9,'(c)','Units','normalized','FontSize',11)
    hold on
    p2 = plot(timed,Rt,'-k');
    plot(timed,ones(length(timed),1), 'color', 'red', 'linewidth',0.3)
    for i = 1:floor(length(indexswitchRt)/2)
        fill([timed(indexswitchRt(1+2*(i-1))) timed(indexswitchRt(2+2*(i-1))) timed(indexswitchRt(2+2*(i-1))) timed(indexswitchRt(1+2*(i-1)))], [0 0 5 5],'r','facealpha',0.075,'EdgeColor','none')
    end
    if length(indexswitchRt)/2 ~= floor(length(indexswitchRt)/2)
            fill([timed(indexswitchRt(end)) timed(end) timed(end) timed(indexswitchRt(end))], [0 0 5 5],'r','facealpha',0.075,'EdgeColor','none')
    end
    ylim([0 5])
%     ylim([min(0.5./Rt) max(2*Rt)])
%     set(gca, 'yscale', 'log')
    yticks([0.2 1 5])
    ylabel('$R_t$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec) 
    set(gca,'Xticklabel',[])
    hold off
    box off
    
    subplot(3,1,3)
    text(0.95,0.9,'(d)','Units','normalized','FontSize',11)
    hold on
    p2 = plot(timed,et,'-k');
        for i = 1:floor(length(indexswitchet)/2)
        fill([timed(indexswitchet(1+2*(i-1))) timed(indexswitchet(2+2*(i-1))) timed(indexswitchet(2+2*(i-1))) timed(indexswitchet(1+2*(i-1)))], [-1 -1 2 2],'r','facealpha',0.075,'EdgeColor','none')
        end
        if length(indexswitchet)/2 ~= floor(length(indexswitchet)/2)
            fill([timed(indexswitchet(end)) timed(end) timed(end) timed(indexswitchet(end))], [-1 -1 2 2],'r','facealpha',0.075,'EdgeColor','none')
        end

        plot(timed,zeros(length(timed),1), 'color', 'red', 'linewidth',0.3)
set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec,'fontsize',11)    
datetick('x','mmm-yy','keeplimits','keepticks')
box off
ylabel('$e_t$ [1/d]', 'fontsize',11, 'interpreter','latex')
%line([timed(1) timed(end)], [0 0], 'color', 'magenta','LineWidth',0.1,'LineStyle','--')
hold off
f.Units='points';
f.Position=[0 0 450 400];
% saveas(f,'plots/norm.pdf','pdf')

