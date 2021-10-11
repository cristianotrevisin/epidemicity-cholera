clear all 
close all
clc

cases_week = csvread('data/cases.csv',1,1)';    %ski the first 2 weeks %from 20 Oct 2010 to 31 Aug 2011 --> the last week is full, so actually to 3 Sep 2011

cases_week = cases_week(:,1:350);

[rSeq1, rSeq2] = load_data();
opt = "best";
[x1,x2] = get50(opt,rSeq1, rSeq2);

[cases_AD1_week0, time] = m4r(x1,x2);

[p1,p2] = get50("95pct",rSeq1, rSeq2);
nopost = 100;
   
ribbon1 = zeros(nopost,350,10);
opt = "rand";
for i = 1:nopost
    i
    [x1,x2] = get50(opt,rSeq1, rSeq2);  
    try
        [cases_AD1_week,time] = m4r(x1,x2);
    catch
        cases_AD1_week = NaN(10,350);
    end

    ribbon1(i,:,:)=cases_AD1_week';
end


    rib_min = zeros(10,350);
    rib_max = zeros(10,350);
    for i = 1:350
        for j = 1:10
            rib_min(j,i) = min(ribbon1(:,i,j));
            rib_max(j,i) = max(ribbon1(:,i,j));
        end
    end


%tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.05.2011','dd.mm.yyyy') ...
%   datenum('01.11.2011','dd.mm.yyyy') datenum('01.05.2012','dd.mm.yyyy') ...
%   datenum('01.11.2012','dd.mm.yyyy') datenum('01.07.2013','dd.mm.yyyy')];

tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.11.2011','dd.mm.yyyy') ...
    datenum('01.11.2012','dd.mm.yyyy') datenum('01.11.2013','dd.mm.yyyy') ...
    datenum('01.11.2014','dd.mm.yyyy') datenum('01.11.2015','dd.mm.yyyy')...
    datenum('01.11.2016','dd.mm.yyyy')];


%%
close all
%close(4000)
ff = figure(4000)
    set(gcf,'color','white')
    ax1 = axes('Position',[0.1 0.4 0.35 0.5]);
    p3 = plot(time,sum(cases_week),'.','color','red');    
    hold on

        p1 = jbfill(time,sum(rib_max),sum(rib_min),[0.410 0.41 0.41],'white','FaceAlpha',0.35);
        hold on
        p2 = plot(time,sum(cases_AD1_week0),'-k');

    box off
    ylabel('Weekly cases')
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    datetick('x','mmm-yy','keeplimits','keepticks')
    set(gca,'Xticklabel',[])
    hold off
    ylim([0 1.5*max(sum(cases_week))])
    title('Haiti', 'fontsize', 20)
    
    ax2 = axes('Position',[0.1 0.1 0.35 0.2]);
    plot(time,cases_week(9,:),'.','color','red');  
    hold on

        jbfill(time,rib_max(9,:),rib_min(9,:),[1 0.647 0],'white','FaceAlpha',0.35);
                hold on
        plot(time,cases_AD1_week0(9,:),'color','#FF8C00')

    box off
    ylim([0 1.5*max(cases_AD1_week0(9,:))])
    %legend('Simulation','Observations', 'fontsize', 14)
    %ylabel('Weekly cases','fontsize',14)
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    datetick('x','mmm-yy','keeplimits','keepticks')
    hold off
    title('Sud')
    
    ax3 = axes('Position',[0.55 0.7 0.35 0.2]);
    plot(time,cases_week(1,:),'.','color','red');   
    hold on
        jbfill(time,rib_max(1,:),rib_min(1,:),[0.117 0.563 1],'white','FaceAlpha',0.35);
        hold on
        plot(time,cases_AD1_week0(1,:),'color','#1E90FF')

    box off
    ylim([0 1.5*max(cases_AD1_week0(1,:))])
    %legend('Simulation','Observations', 'fontsize', 14)
    %ylabel('Weekly cases','fontsize',14)
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    datetick('x','mmm-yy','keeplimits','keepticks')
    set(gca,'Xticklabel',[])
    hold off
    title('Artibonite')
    
    ax4 = axes('Position',[0.55 0.4 0.35 0.2]);
    plot(time,cases_week(2,:),'.','color','red');
    hold on

        jbfill(time,rib_max(2,:),rib_min(2,:),[1 0.410 0.703],'white','FaceAlpha',0.35);
                hold on
        plot(time,cases_AD1_week0(2,:),'color','#FF69B4')

    box off
    ylim([0 1.5*max(cases_AD1_week0(2,:))])
    %legend('Simulation','Observations', 'fontsize', 14)
    %ylabel('Weekly cases','fontsize',14)
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    datetick('x','mmm-yy','keeplimits','keepticks')
    set(gca,'Xticklabel',[])
    hold off
    title('Centre')
    
    ax5 = axes('Position',[0.55 0.1 0.35 0.2]);
    plot(time,cases_week(8,:),'.','color','red');
    hold on

        jbfill(time,rib_max(8,:),rib_min(8,:),[0 1 0],'white','FaceAlpha',0.35);
                hold on
        plot(time,cases_AD1_week0(8,:),'color','#00FF00')

    box off
    ylim([0 1.5*max(cases_AD1_week0(8,:))])
    %legend('Simulation','Observations', 'fontsize', 14)
    %ylabel('Weekly cases','fontsize',14)
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    datetick('x','mmm-yy','keeplimits','keepticks')
    hold off
    title('Ouest')
ff.Units='points';
ff.Position=[0 0 450 400];  
saveas(ff,'fitting.pdf','pdf')
%     ax6 = axes('Position',[.17 .55 .30 .30]);
% box on;
% fileName = 'maps4plot.png';
% rgbImage = imread(fileName);
% imshow(rgbImage);
% grid off;