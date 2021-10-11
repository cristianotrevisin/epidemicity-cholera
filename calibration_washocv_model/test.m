clc
close all
clear all

cases_week = csvread('data/cases.csv',1,1)';    %ski the first 2 weeks %from 20 Oct 2010 to 31 Aug 2011 --> the last week is full, so actually to 3 Sep 2011

cases_week = cases_week(:,1:350); 

vec = [1 3 0.1 0.1 2 2];
[ModPred,times_data] = m5c(vec);

tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.11.2011','dd.mm.yyyy') ...
    datenum('01.11.2012','dd.mm.yyyy') datenum('01.11.2013','dd.mm.yyyy') ...
    datenum('01.11.2014','dd.mm.yyyy') datenum('01.11.2015','dd.mm.yyyy')...
    datenum('01.11.2016','dd.mm.yyyy') datenum('01.11.2017','dd.mm.yyyy')...
    datenum('01.11.2018','dd.mm.yyyy')];

figure(4002)
    set(gcf,'color','white')
    p3 = plot(times_data,sum(cases_week(:,1:end)),'o','markeredgecolor','red','markerfacecolor','white','markersize',4.5,'color',[0.5 0.5 0.5],'linewidth', 1);    
    hold on
    p2 = plot(times_data,sum(ModPred),'linewidth',1.5,'color','blue');
    box off
    legend([p2 p3], 'Model','Observations', 'fontsize', 14)
    ylabel('Weekly cases','fontsize',14)
    set(gca, 'fontsize', 14)
    set(gca,'Xlim',[times_data(1) times_data(end)],'Xtick',tick_vec)
    datetick('x','mmm-yy','keeplimits','keepticks')
    hold off