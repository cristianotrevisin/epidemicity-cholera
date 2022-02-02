%%% COMPARISON BETWEEN RT AND ET (RATIO)
% 

close all
clear all
clc

[rSeq1, rSeq2] = load_data();
opt = "best";
[x1,x2] = get50(opt,rSeq1, rSeq2);
%%
[cases_AD1_week, time, y, cati, ocv] = SIARBV(2, 2, x1,x2);
[Rta, eta,ratio] = diagnosis(y, 2, 2, x1, x2);

et = smoothdata(eta,'movmean',28);
Rt = smoothdata(Rta,'movmean',28);


%% PLOT RATIO
timed = datenum('2010-10-20'):1:datenum('2017-07-01');
tick_vec=[datenum('01.11.2010','dd.mm.yyyy') datenum('01.11.2011','dd.mm.yyyy') ...
    datenum('01.11.2012','dd.mm.yyyy') datenum('01.11.2013','dd.mm.yyyy') ...
    datenum('01.11.2014','dd.mm.yyyy') datenum('01.11.2015','dd.mm.yyyy')...
    datenum('01.11.2016','dd.mm.yyyy')];

fb = figure(9999)
    hold on
    plot(timed,1./ratio,'Color', [0 0.4470 0.7410],'Linewidth',0.005)
    plot(timed,ones(length(timed),1)*mean(1./ratio), 'color','red')
    set(gca,'Xlim',[timed(1) timed(end)],'Xtick',tick_vec,'fontsize',11)    
    datetick('x','mmm-yy','keeplimits','keepticks')
    box off
    ylabel('$\zeta_t/\epsilon_t = \mathcal{R}_{t,lim}$', 'fontsize',11, 'interpreter','latex')
    legend('Values','Mean')
    ylim([0.72 0.74])
    hold off
    fb.Units='points';
    fb.Position=[0 0 450 200];

%% PRE-SMOOTHING
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
Rsq1 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)

%% SMOOTHED
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