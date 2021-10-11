clc
close all
clear all

% load data
load ../markov_chains/TEMP_RES_2.mat Sequences output

% remove zeros
index0=find(squeeze(Sequences(:,end-1,1)),1,'first');
index1=find(squeeze(Sequences(:,1,1)),1,'last');
Sequences =  Sequences(index0:index1,:,:);

% INITIALISE
tailleng = 1000000*3;
for i = 1:size(Sequences,2)
        Seq(:,i) = reshape([Sequences(:,i,1)'; Sequences(:,i,2)'; Sequences(:,i,3)'], [], 1)';
   end

MCMCPar.n = size(Sequences,2)-2;
titlestr = ["\xi_1", "\xi_2", "t_1", "t_2"];

% MARKOV CHAIN PLOT
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

% CONVERGENCE DIAGNOSTIC
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



% POSTERIOR DISTRIBUTION
figure(5555)
    for i=1:MCMCPar.n
        subplot(2,ceil((MCMCPar.n)/2),i)
        histogram(Seq(end-tailleng:end,i))
        %set(gca,'xlim',[ParRange.minn(i) ParRange.maxn(i)])
        set(gca,'fontsize',14)
        title(titlestr(i));
    end

