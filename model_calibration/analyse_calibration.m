clc
close all
clear all

% load data
load ../markov_chains/TEMP_RES_1.mat Sequences output

% remove zeros
index0=find(squeeze(Sequences(:,end-1,1)),1,'first');
index1=find(squeeze(Sequences(:,1,1)),1,'last');
Sequences1 =  Sequences(index0:index1,:,:);
index2=find(output.R_stat(:,1),1,'last');
R_stat1 = output.R_stat(1:index2,:);
MCMCPar.n1 = size(Sequences,2)-2;

clear Sequences output index0 index1 index2

% load data
load ../markov_chains/TEMP_RES_2.mat Sequences output

% remove zeros
index0=find(squeeze(Sequences(:,end-1,1)),1,'first');
index1=find(squeeze(Sequences(:,1,1)),1,'last');
Sequences2 =  Sequences(index0:index1,:,:);
index2=find(output.R_stat(:,1),1,'last');
R_stat2 = output.R_stat(1:index2,:);
MCMCPar.n2 = size(Sequences,2)-2;

clear Sequences output index0 index1 index2

titlestr1 = ["\theta", "m", "D", "\phi", "\rho", "\sigma", "\mu_B", "\beta_0", "\psi", "t_0", "r"];
titlestr2 = ["\xi_1", "\xi_2", "t_w"];

%% MARKOV CHAIN PLOT
f = figure(2222)
    for i=1:MCMCPar.n1
        subplot(5,3,i)
        plot(squeeze(Sequences1(:,i,:)))
        set(gca,'xlim',[0 size(Sequences1,1)*1.2])
        set(gca,'fontsize',11)
        title(titlestr1(i));
    end
    for j=1:MCMCPar.n2
        subplot(5,3,i+j)
        plot(squeeze(Sequences2(:,j,:)))
        set(gca,'xlim',[0 size(Sequences2,1)*1.2])
        set(gca,'fontsize',11)
        title(titlestr2(j));
    end
f.Units='points';
f.Position=[0 0 400 600];
    
    

%% CONVERGENCE DIAGNOSTIC
g = figure(3333)
subplot(2,1,1)
hold on
for i=1:MCMCPar.n1
    plot(R_stat1(:,1), R_stat1(:,i+1))
end
ylabel('$\hat{R}$','Interpreter', 'latex','fontsize',11)
set(gca,'yscale','log')
set(gca,'YMinorTick','on')
line([1 R_stat1(end,1)], [1.05 1.05], 'Color','red','LineStyle','--','LineWidth',1)
legend([titlestr1])
legend boxoff
ylim([1 10])
grid on

subplot(2,1,2)
hold on
for i=1:MCMCPar.n2
    plot(R_stat2(:,1), R_stat2(:,i+1))
end
set(gca,'yscale','log')
set(gca,'YMinorTick','on')
line([1 R_stat2(end,1)], [1.05 1.05], 'Color','red','LineStyle','--','LineWidth',1)
ylabel('$\hat{R}$','Interpreter', 'latex','fontsize',11)
legend([titlestr2])
legend boxoff
ylim([1 10])
grid on
xlabel('Iteration','fontsize',11)
g.Units='points';
g.Position=[0 0 400 400];

%% POSTERIOR DISTRIBUTIONS

[rSeq1, rSeq2] = load_data();
[xb1,xb2] = get50("best",rSeq1, rSeq2);
[x51,x52] = get50("95pct",rSeq1, rSeq2);
% POSTERIOR DISTRIBUTION
h = figure(5555)
    for i=1:MCMCPar.n1
        subplot(5,3,i)
        histogram(rSeq1(:,i),'EdgeColor','none')
        set(gca,'fontsize',11)
        xline(xb1(i),'--k','Linewidth',1)
        xline(x51(i,1),'--r','Linewidth',0.5); xline(x51(i,2),'--r','Linewidth',0.5)
        title(titlestr1(i));
    end
    for j=1:MCMCPar.n2
        subplot(5,3,i+j)
        histogram(rSeq2(:,j),'EdgeColor','none')
        set(gca,'fontsize',11)
        xline(xb2(j),'--k','Linewidth',1)
        xline(x52(j,1),'--r','Linewidth',0.5); xline(x52(j,2),'--r','Linewidth',0.5)
        title(titlestr2(j));
    end
h.Units='points';
h.Position=[0 0 400 600];