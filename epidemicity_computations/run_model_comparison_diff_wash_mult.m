close all
clear all
clc

load data/geodata POPnodes WS_dept
POPnodes = POPnodes' * WS_dept;
POPnodes = POPnodes';



cases_week = csvread('data/cases.csv',1,1)';
cases_week = cases_week(:,1:350); 
multipliers = [0 0.5 1 2 5 10];
cases = zeros(length(multipliers),350);
Rt = zeros(length(multipliers),2447);
et = zeros(length(multipliers),2447);
for i = 1:length(multipliers)
    [cases_AD1_week, time, y] = m5c_washmult(multipliers(i));
    cases(i,:) = sum(cases_AD1_week,1);
    [a, b] = diagnosis(y,multipliers(i),1);
    Rt(i,:) = smoothdata(a,'movmean',28);
    et(i,:) = smoothdata(b,'movmean',28);
end



%% PLOT
close all
timed = datenum('2010-10-20'):1:datenum('2017-07-01');
tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.11.2011','dd.mm.yyyy') ...
    datenum('01.11.2012','dd.mm.yyyy') datenum('01.11.2013','dd.mm.yyyy') ...
    datenum('01.11.2014','dd.mm.yyyy') datenum('01.11.2015','dd.mm.yyyy')...
    datenum('01.11.2016','dd.mm.yyyy')];

tickvec = [-0.05, 0, 1, 100];

%%

f = figure(4002);

    subplot(3,1,1) 
    text(0.95,0.9,'(a)','Units','normalized','FontSize',11)
    hold on
    set(gcf,'color','white')
    p2 = bar(time, sum(cases_week),'FaceColor','#C1C1C1')
    for i = 1:length(multipliers)
        plot(time,cases(i,:));
    end
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
    text(0.95,0.9,'(b)','Units','normalized','FontSize',11)
    hold on
    for i = 1:length(multipliers)
    plot(timed,Rt(i,:));
    end
    lgd = legend('0', '0.5', '1', '2', '5', '10','NumColumns',2);
    title(lgd, 'WaSH multiplier')
    h=plot(timed,ones(length(timed),1), 'color', 'red', 'linewidth',0.3)
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';

%     ylim([min(0.5./Rt) max(2*Rt)])
%     set(gca, 'yscale', 'log')
    yticks([0.2 1 5])
    ylabel('$R_t$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec) 
    set(gca,'Xticklabel',[])
    hold off
    box off
    
    subplot(3,1,3)
    text(0.95,0.9,'(c)','Units','normalized','FontSize',11)
    hold on
    for i = 1:length(multipliers)
        plot(timed,et(i,:));
    end
plot(timed,zeros(length(timed),1), 'color', 'red', 'linewidth',0.3)

set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec,'fontsize',11)    
datetick('x','mmm-yy','keeplimits','keepticks')
box off
%line([timed(1) timed(end)], [0 0], 'color', 'magenta','LineWidth',0.1,'LineStyle','--')
hold off
f.Units='points';
f.Position=[0 0 450 400];
% saveas(f,'plots/norm.pdf','pdf')

