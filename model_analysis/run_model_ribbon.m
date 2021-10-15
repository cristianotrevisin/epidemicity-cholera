close all
clear all
clc

load ../data/geodata POPnodes WS_dept
POPnodes = POPnodes' * WS_dept;
POPnodes = POPnodes';

[rSeq1, rSeq2] = load_data();
opt = "best";
[x1,x2] = get50(opt,rSeq1, rSeq2);

cases_week = csvread('../data/cases.csv',1,1)';
cases_week = cases_week(:,1:350); 

[cases_AD1_week0, time, y, cati_sum, ocv, rain] = SIARBV(2, 2. x1,x2);
[Rt0, et0] = diagnosis(y, 2, 2, x1, x2);
NS = 1 - sum((cases_AD1_week0(:) - cases_week(:)).^2)/sum((cases_week(:) - mean(cases_week(:))).^2)
et0 = smoothdata(et0,'movmean',28);
Rt0 = smoothdata(Rt0,'movmean',28);


nopost = 5;
   
ribbon1cases = zeros(nopost,350);
ribbon1norm = zeros(nopost,length(et0));
ribbon1et = zeros(nopost,length(et0));
ribbon1Rt = zeros(nopost,length(et0));


opt = "rand";
for j = 1:nopost
    j
    [x1,x2] = get50(opt,rSeq1, rSeq2);
%     [x3,x4] = get90(opt,rSeq1, rSeq2);

    try
        [cases_AD1_week, ~, y, ~, ] = m4r(x1,x2);
        [Rt, et] = diagnosis4r(y, x1, x2);
                
        et = smoothdata(et,'movmean',28);
        Rt = smoothdata(Rt,'movmean',28);
    catch
        
    end

    
    sumc = sum(cases_AD1_week,1);
    ribbon1cases(j,:)= sumc;
    ribbon1et(j,:) = et;
    ribbon1Rt(j,:) = Rt;
end


    rib_min_c = zeros(1,350);
    rib_max_c = zeros(1,350);
    rib_min_e = zeros(1,length(et0));
    rib_max_e = zeros(1,length(et0));
    rib_min_r = zeros(1,length(et0));
    rib_max_r = zeros(1,length(et0));
    for i = 1:350
        rib_min_c(i) = min(ribbon1cases(:,i));
        rib_max_c(i) = max(ribbon1cases(:,i));
    end
    for i = 1:length(et0)
        rib_min_e(i) = min(ribbon1et(:,i));
        rib_max_e(i) = max(ribbon1et(:,i));
        rib_min_r(i) = min(ribbon1Rt(:,i));
        rib_max_r(i) = max(ribbon1Rt(:,i));
    end

%% PLOT
close all
timed = datenum('2010-10-20'):1:datenum('2017-07-01');
tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.11.2011','dd.mm.yyyy') ...
    datenum('01.11.2012','dd.mm.yyyy') datenum('01.11.2013','dd.mm.yyyy') ...
    datenum('01.11.2014','dd.mm.yyyy') datenum('01.11.2015','dd.mm.yyyy')...
    datenum('01.11.2016','dd.mm.yyyy')];

tick_vec1=[datenum('01.01.2011','dd.mm.yyyy') datenum('01.01.2012','dd.mm.yyyy') ...
    datenum('01.01.2013','dd.mm.yyyy') datenum('01.01.2014','dd.mm.yyyy') ...
    datenum('01.01.2015','dd.mm.yyyy') datenum('01.01.2016','dd.mm.yyyy')...
    datenum('01.01.2017','dd.mm.yyyy')];

indexswitchRt = [1 find(diff(sign(Rt0-1)))+1];
indexswitchet = [1 find(diff(sign(et0)))+1];

tickvec = [-0.05, 0, 1, 100];
f = figure(4002);

    ax1 = axes('Position',[0.1 0.55 0.8 0.35]);
    text(0.95,0.2,'(a)','Units','normalized','FontSize',11)
    hold on
    set(gcf,'color','white')
    plot(time,sum(cases_week),'.','color','red');    

    jbfill(time,rib_max_c,rib_min_c,[0.410 0.41 0.41],'white','FaceAlpha',0.5);
        hold on
        p2 = plot(time,sum(cases_AD1_week0),'-k');
    box off
    %p2 = plot(time,sum(cases_week),'o','markeredgecolor','red','markerfacecolor','white','markersize',2,'color',[0.5 0.5 0.5],'linewidth', 1);
    %legend('Observation', 'Simulation','fontsize', 11, 'location','southeast')
    ylabel('$\Delta C$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    set(gca,'Xticklabel',[])
    %legend boxoff
    box off
    ylim([0 30000])
    hold off

%     ax2 = axes('Position',[0.35 0.75 0.5 0.15]);
%     text(0.95,0.9,'(b)','Units','normalized','FontSize',11)
%     hold on
%     jbfill(timed,rib_max_i,rib_min_i,[0.410 0.41 0.41],'white','FaceAlpha',0.35);
%         hold on
%         p2 = plot(timed,i0/i0(1),'-k');
%     ylabel('$\frac{||\bf{y}(t)||}{||\bf{y}(0)||}$', 'fontsize',11, 'Interpreter','latex')
%     set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec)   
%     set(gca,'Xticklabel',[])
%     ylim([0 20])
%     hold off
%     box off
    cati_d = diff([zeros(1,10) ; cati_sum]);
    ax2 = axes('Position',[0.315 0.75 0.5 0.15]);
    text(0.95,0.9,'(b)','Units','normalized','FontSize',11)
    yyaxis left
    hold on
    bar(timed, cati_d, 'stacked')
    ylabel('Weekly NPI')
    set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec1)   
    datetick('x','yyyy','keeplimits','keepticks')
    yyaxis right
    set(gca,'Ydir','reverse')
    bar(timed, rain, 'b')
    ylabel('$J$ [mm/d]', 'fontsize',11, 'interpreter','latex')
    ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
    hold off
    box off

    ax3 = axes('Position',[0.1 0.325 0.8 0.175]);
    text(0.95,0.9,'(c)','Units','normalized','FontSize',11)
    hold on
    set(gca,'Xlim',[timed(1) timed(end)])   
    jbfill(timed,rib_max_r,rib_min_r,[0.410 0.41 0.41],'white','FaceAlpha',0.35);
        hold on
    p2 = plot(timed,Rt0,'-k');
    plot(timed,ones(length(timed),1), 'color', 'red', 'linewidth',0.3)
    for i = 1:floor(length(indexswitchRt)/2)
        fill([timed(indexswitchRt(1+2*(i-1))) timed(indexswitchRt(2+2*(i-1))) timed(indexswitchRt(2+2*(i-1))) timed(indexswitchRt(1+2*(i-1)))], [0 0 5 5],'r','facealpha',0.075,'EdgeColor','none')
    end
    if length(indexswitchRt)/2 ~= floor(length(indexswitchRt)/2)
            fill([timed(indexswitchRt(end)) timed(end) timed(end) timed(indexswitchRt(end))], [0 0 5 5],'r','facealpha',0.075,'EdgeColor','none')
    end
    ylim([0.1 3])
    ylabel('$\mathcal{R}_t$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec) 
    set(gca,'Xticklabel',[])
    hold off
    box off
    
    ax4 = axes('Position',[0.1 0.1 0.8 0.175]);

    text(0.95,0.9,'(d)','Units','normalized','FontSize',11)
    hold on
    set(gca,'Xlim',[timed(1) timed(end)])
    plot(timed,zeros(length(timed),1), 'color', 'red', 'linewidth',0.3)
    jbfill(timed,rib_max_e,rib_min_e,[0.410 0.41 0.41],'white','FaceAlpha',0.35);
        hold on
        p2 = plot(timed,et0,'-k');
        for i = 1:floor(length(indexswitchet)/2)
        fill([timed(indexswitchet(1+2*(i-1))) timed(indexswitchet(2+2*(i-1))) timed(indexswitchet(2+2*(i-1))) timed(indexswitchet(1+2*(i-1)))], [-1 -1 5 5],'r','facealpha',0.075,'EdgeColor','none')
        end
        if length(indexswitchet)/2 ~= floor(length(indexswitchet)/2)
            fill([timed(indexswitchet(end)) timed(end) timed(end) timed(indexswitchet(end))], [-1 -1 5 5],'r','facealpha',0.075,'EdgeColor','none')
        end
        ylabel('$e_t$ [1/d]', 'fontsize',11, 'interpreter','latex')     
    ylim([-0.2 1])
set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec,'fontsize',11)    
datetick('x','mmm-yy','keeplimits','keepticks')
box off
%line([timed(1) timed(end)], [0 0], 'color', 'magenta','LineWidth',0.1,'LineStyle','--')
hold off
f.Units='points';
f.Position=[0 0 400 400];

%saveas(f,'plots/normm.pdf','pdf')