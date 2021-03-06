clear all 
close all
clc

[rSeq1, rSeq2] = load_data();
opt = "best";
[x1,x2] = get50(opt,rSeq1, rSeq2);

[cases_AD1_week, time] = SIARBV(2, 2, x1,x2);
cases_AD1_week = cases_AD1_week(:,141:end);
baseline = sum(sum(cases_AD1_week));

for i = 1:11
    x1(i) = x1(i)*1.1;
    [cases_AD1_week] = SIARBV(2, 2, x1,x2);
    cases_AD1_week = cases_AD1_week(:,141:end);
    up(i) = sum(sum(cases_AD1_week))/baseline-1;
    x1(i) = x1(i)/1.1/1.1;
    [cases_AD1_week] = SIARBV(2, 2, x1,x2);
    cases_AD1_week = cases_AD1_week(:,141:end);
    down(i) = sum(sum(cases_AD1_week))/baseline-1;
    x1(i) = x1(i)*1.1;
end

for i = 1:3
    x2(i) = x2(i)*1.1;
    [cases_AD1_week, time] = SIARBV(2, 2, x1,x2);
    cases_AD1_week = cases_AD1_week(:,141:end);
    up(11+i) = sum(sum(cases_AD1_week))/baseline-1;
    x2(i) = x2(i)/1.1/1.1;
    [cases_AD1_week, time] = SIARBV(2, 2, x1,x2);
    cases_AD1_week = cases_AD1_week(:,141:end);
    down(11+i) = sum(sum(cases_AD1_week))/baseline-1;
    x2(i) = x2(i)*1.1;
end


%% plot
close all
f = figure(92648)
barh(1:14,down*100)
hold on
barh(1:14,up*100)
xlabel('[%] prevalence variation')
box off
set(gca,'Ydir','reverse')
ax = gca;
%set(h,'ycolor','w')
ax.YColor = 'w';
ax.YAxis.Label.Color='k';
titlestr = ["\theta", "m", "D", "\phi", "\rho", "\sigma", "\mu_B", "\beta_0", "\psi", "t_0", "r", "\xi_1", "\xi_2", "t_w"];
yticks(1:15)
yticklabels(titlestr)
f.Units='points';
f.Position=[0 0 400 350];
saveas(f,'sensitivity.pdf','pdf')