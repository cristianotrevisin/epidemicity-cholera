clear all
close all
clc
load '../data/ocv.mat' datespace rv_1d rv_2d eta_1d eta_2d
eta_1d(:,4:7)=[]; eta_2d(:,4:7)=[]; rv_1d(:,4:7)=[]; rv_2d(:,4:7)=[];
eta_1d(:,6)=[]; eta_2d(:,6)=[]; rv_1d(:,6)=[]; rv_2d(:,6)=[];
a = datenum(2016,07,01);
b = datenum(2019,01,31);

f=figure()
subplot(2,2,1)
    hold on
    for i = 1:5
        stairs(datespace,cumsum(rv_1d(:,i)),'linewidth',1)
    end
    set(gca,'Xlim',[a b])
    datetick('x','mmm-yy','keeplimits','keepticks')
    box off
    title('First doses')
    ylabel('$v^d + w^d$', 'interpreter','latex')
subplot(2,2,2)
    hold on
    for i = 1:5
        stairs(datespace,cumsum(rv_2d(:,i)),'linewidth',1)
    end
    set(gca,'Xlim',[a b])
    datetick('x','mmm-yy','keeplimits','keepticks')
    box off
    title('Second doses')
    ylabel('$v^{dd}$', 'interpreter','latex')
subplot(2,2,3)
    hold on
    for i = 1:5
        plot(datespace,eta_1d(:,i),'linewidth',1)
    end
    set(gca,'Xlim',[a b])
    datetick('x','mmm-yy','keeplimits','keepticks')
    box off
    title('Efficacy 1st dose')
    xlabel('Time')
    ylabel('$\eta^d$', 'interpreter','latex')
subplot(2,2,4)
    hold on
    for i = 1:5
        plot(datespace,eta_2d(:,i),'linewidth',1)
    end
    set(gca,'Xlim',[a b])
    title('Efficacy 2nd dose')
    ylabel('$\eta^{dd}$', 'interpreter','latex')
    xlabel('Time')
    datetick('x','mmm-yy','keeplimits','keepticks')
    legend("Artibonite", "Centre", "Grande Anse", "Ouest", "Sud", "Location","bestoutside")
    legend boxoff
    box off
    f.Units='points';
    f.Position=[0 0 400 400];