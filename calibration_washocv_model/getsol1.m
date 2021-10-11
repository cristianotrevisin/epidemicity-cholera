clear all
close all
clc

load TEMP_RES.mat Sequences output

% close all
index0=find(squeeze(Sequences(:,end-1,1)),1,'first');
index1=find(squeeze(Sequences(:,1,1)),1,'last');

ribbon = 0;
Sequences =  Sequences(index0:index1,:,:);

[iter, chain] = find(squeeze(Sequences(:,end-1,:))==max(max(squeeze(Sequences(:,end-1,:)))));

vec = Sequences(max(iter),1:end-2,max(chain))

cases_week = csvread('data/cases.csv',1,1)';    %ski the first 2 weeks %from 20 Oct 2010 to 31 Aug 2011 --> the last week is full, so actually to 3 Sep 2011

cases_week = cases_week(:,1:350);

ModelName = 'm5c';

evalstr = ['[cases_AD1_week, time] = ',ModelName,'(vec);']; eval(evalstr)

NS = 1 - sum((cases_AD1_week(:) - cases_week(:)).^2)/sum((cases_week(:) - mean(cases_week(:))).^2)
tailleng = 3*25000;
if ribbon
   nopost = 1000;
   ribbon1 = zeros(nopost,141,10);
   clear Seq
   for i = 1:size(Sequences,2)
        Seq(:,i) = reshape([Sequences(:,i,1)'; Sequences(:,i,2)'; Sequences(:,i,3)'], [], 1)';
   end
   
   for i = 1:nopost
       i
       k = randi([1 tailleng]);
       vec = Seq(end-k,:);
       evalstr = ['[cases_AD1_week,time] = ',ModelName,'(vec);']; eval(evalstr)
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
end

%tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.05.2011','dd.mm.yyyy') ...
%   datenum('01.11.2011','dd.mm.yyyy') datenum('01.05.2012','dd.mm.yyyy') ...
%   datenum('01.11.2012','dd.mm.yyyy') datenum('01.07.2013','dd.mm.yyyy')];

tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.11.2011','dd.mm.yyyy') ...
    datenum('01.11.2012','dd.mm.yyyy') datenum('01.11.2013','dd.mm.yyyy') ...
    datenum('01.11.2014','dd.mm.yyyy') datenum('01.11.2015','dd.mm.yyyy')...
    datenum('01.11.2016','dd.mm.yyyy')];
%%
%close(4000)
figure(4000)
    set(gcf,'color','white')
    ax1 = axes('Position',[0.1 0.4 0.35 0.5]);
    p3 = plot(time,sum(cases_week),'o','markeredgecolor','red','markerfacecolor','white','markersize',4.5,'color',[0.5 0.5 0.5],'linewidth', 1);    
    hold on
    if ribbon==1
        p1 = jbfill(time,sum(rib_max),sum(rib_min),[0 0.5 0.8],'white','FaceAlpha',0.35);
        hold on
        p2 = plot(time,sum(cases_AD1_week),'-b');
    else
        p2 = plot(time,sum(cases_AD1_week),'-b','linewidth',1.5);
    end
    box off
    legend([p2 p3], 'Simulation','Observations', 'fontsize', 14)
    ylabel('Weekly cases','fontsize',14)
    set(gca, 'fontsize', 14)
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    datetick('x','mmm-yy','keeplimits','keepticks')
    %set(gca,'Xticklabel',[])
    hold off
    title('Haiti', 'fontsize', 20)
    
    ax2 = axes('Position',[0.1 0.1 0.35 0.2]);
    plot(time,cases_week(9,:),'o','markeredgecolor','red','markerfacecolor','white','markersize',4.5,'color',[0.5 0.5 0.5],'linewidth', 1);  
    hold on
    if ribbon
        jbfill(time,rib_max(9,:),rib_min(9,:),[0 0.5 0.8],'white','FaceAlpha',0.35);
                hold on
        plot(time,cases_AD1_week(9,:),'-b')
    else
        plot(time,cases_AD1_week(9,:),'-b','linewidth',1.5)
        
    end
    box off
    %legend('Simulation','Observations', 'fontsize', 14)
    %ylabel('Weekly cases','fontsize',14)
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    set(gca, 'fontsize', 14)
    datetick('x','mmm-yy','keeplimits','keepticks')
    hold off
    title('Sud')
    
    ax3 = axes('Position',[0.55 0.7 0.35 0.2]);
    plot(time,cases_week(1,:),'o','markeredgecolor','red','markerfacecolor','white','markersize',4.5,'color',[0.5 0.5 0.5],'linewidth', 1);   
    hold on
    if ribbon
        jbfill(time,rib_max(1,:),rib_min(1,:),[0 0.5 0.8],'white','FaceAlpha',0.35);
        hold on
        plot(time,cases_AD1_week(1,:),'-b')
    else
        plot(time,cases_AD1_week(1,:),'-b','linewidth',1.5)
        
    end
    box off
    %legend('Simulation','Observations', 'fontsize', 14)
    %ylabel('Weekly cases','fontsize',14)
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    set(gca, 'fontsize', 14)
    datetick('x','mmm-yy','keeplimits','keepticks')
    hold off
    title('Artibonite')
    
    ax4 = axes('Position',[0.55 0.4 0.35 0.2]);
    plot(time,cases_week(2,:),'o','markeredgecolor','red','markerfacecolor','white','markersize',4.5,'color','white','markersize',4.5,'color',[0.5 0.5 0.5],'linewidth', 1);
    hold on
    if ribbon
        jbfill(time,rib_max(2,:),rib_min(2,:),[0 0.5 0.8],'white','FaceAlpha',0.35);
                hold on
        plot(time,cases_AD1_week(2,:),'-b')
    else
        plot(time,cases_AD1_week(2,:),'-b','linewidth',1.5)
        
    end
    box off
    %legend('Simulation','Observations', 'fontsize', 14)
    %ylabel('Weekly cases','fontsize',14)
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    set(gca, 'fontsize', 14)
    datetick('x','mmm-yy','keeplimits','keepticks')
    hold off
    title('Centre')
    
    ax5 = axes('Position',[0.55 0.1 0.35 0.2]);
    plot(time,cases_week(8,:),'o','markeredgecolor','red','markerfacecolor','white','markersize',4.5,'color',[0.5 0.5 0.5],'linewidth', 1);
    hold on
    if ribbon
        jbfill(time,rib_max(8,:),rib_min(8,:),[0 0.5 0.8],'white','FaceAlpha',0.35);
                hold on
        plot(time,cases_AD1_week(8,:),'-b')
    else
        plot(time,cases_AD1_week(8,:),'-b','linewidth',1.5)
       
    end
    box off
    %legend('Simulation','Observations', 'fontsize', 14)
    %ylabel('Weekly cases','fontsize',14)
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    set(gca, 'fontsize', 14)
    datetick('x','mmm-yy','keeplimits','keepticks')
    hold off
    title('Ouest')
    
    ax6 = axes('Position',[.22 .58 .20 .20]);
box on;
fileName = 'map.png';
rgbImage = imread(fileName);
imshow(rgbImage);
grid off;
%%    
edit =  1  ;

if edit
    MCMCPar.n = size(Sequences,2)-2;
    titlestr = ["b_1", "b_2", "t_1", "t_2"];

% titlestr = ["\lambda_1", "\lambda_2", "t_1"];
    
    figure(2222)


    for i=1:MCMCPar.n
        subplot(3,ceil((MCMCPar.n+1)/3),i)
        plot(squeeze(Sequences(:,i,:)))
        set(gca,'xlim',[0 size(Sequences,1)*1.2])
        set(gca,'fontsize',14)
        title(titlestr(i));
    end
    subplot(3,ceil((MCMCPar.n+1)/3),i+1)
    semilogy(squeeze(Sequences(:,i+1,:)))
    title('like')
    
    index2=find(output.R_stat(:,1),1,'last');

    figure(3333)
    hold on
    for i=1:MCMCPar.n
        plot(output.R_stat(1:index2,1), output.R_stat(1:index2,i+1))
    end
    ylabel('$\hat{R}$','Interpreter', 'latex')
    set(gca,'yscale','log')
    set(gca,'YMinorTick','on')
    line([1 output.R_stat(index2,1)], [1.05 1.05], 'Color','red','LineStyle','--','LineWidth',2)
    legend([titlestr])
    xlabel('Iter')
    legend boxoff
    ylim([1 10])
end



%%
posterior = 1;
if posterior
     for i = 1:size(Sequences,2)
        Seq(:,i) = reshape([Sequences(:,i,1)'; Sequences(:,i,2)'; Sequences(:,i,3)'], [], 1)';
   end

    figure(5555)

    for i=1:MCMCPar.n
        subplot(2,ceil((MCMCPar.n)/2),i)
        histogram(Seq(end-tailleng:end,i))
        %set(gca,'xlim',[ParRange.minn(i) ParRange.maxn(i)])
        set(gca,'fontsize',14)
        title(titlestr(i));
    end
end
