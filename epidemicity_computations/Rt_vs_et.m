close all
clear all
clc

load data/geodata POPnodes WS_dept
POPnodes = POPnodes' * WS_dept;
POPnodes = POPnodes';

[rSeq1, rSeq2] = load_data();
opt = "best";
[x1,x2] = get50(opt,rSeq1, rSeq2);

cases_week = csvread('data/cases.csv',1,1)';
cases_week = cases_week(:,1:350); 

[cases_AD1_week, time, y, cati_sum, ocv, rain] = m4r(x1,x2);
[Rta, eta, rall, eall] = diagnosis4r(y, x1, x2);

et = smoothdata(eta,'movmean',28);
Rt = smoothdata(Rta,'movmean',28);


%% PLOT
close all
timed = datenum('2010-10-20'):1:datenum('2017-07-01');
tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.11.2011','dd.mm.yyyy') ...
    datenum('01.11.2012','dd.mm.yyyy') datenum('01.11.2013','dd.mm.yyyy') ...
    datenum('01.11.2014','dd.mm.yyyy') datenum('01.11.2015','dd.mm.yyyy')...
    datenum('01.11.2016','dd.mm.yyyy')];

tickvec = [-0.05, 0, 1, 100];



f = figure(4002);

    subplot(3,1,1) 
    text(0.95,0.9,'(a)','Units','normalized','FontSize',11)
    hold on
    set(gcf,'color','white')
    p2 = bar(time, sum(cases_week),'FaceColor','#C1C1C1')
    plot(time,sum(cases_AD1_week),'color','black','Linewidth',1);
    box off
    ylabel('$\Delta C$', 'fontsize',11, 'interpreter','latex')
    set(gca,'Xlim',[time(1) time(end)],'Xtick',tick_vec)
    set(gca,'Xticklabel',[])
    box off
    ylim([0 26000])
    hold off

    subplot(3,1,2)
    text(0.95,0.9,'(c)','Units','normalized','FontSize',11)
    hold on
    p2 = plot(timed,Rt,'-k');
    plot(timed,ones(length(timed),1), 'color', 'red', 'linewidth',0.3)
        ylim([0 5])
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
plot(timed,zeros(length(timed),1), 'color', 'red', 'linewidth',0.3)
set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec,'fontsize',11)    
datetick('x','mmm-yy','keeplimits','keepticks')
box off
ylabel('$e_t$ [1/d]', 'fontsize',11, 'interpreter','latex')
hold off
f.Units='points';
f.Position=[0 0 450 400];


%%
x = Rta';
y = eta';


X = [ones(length(x),1) x];
b = X\y

yCalc2 = X*b;


figure()
scatter(x,y)
hold on
xlabel('Rt')
ylabel('et')
grid on
plot(x,yCalc2,'--')
legend('Data','Fit','Location','best');
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)

%%
x = Rt';
y = et';


X = [ones(length(x),1) x];
b = X\y

yCalc2 = X*b;


figure()
scatter(x,y)
hold on
xlabel('Rt')
ylabel('et')
grid on
plot(x,yCalc2,'--')
legend('Data','Fit','Location','best');
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)

%%
figure()
scatter(rall, eall)