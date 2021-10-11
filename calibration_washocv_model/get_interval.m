%function [x] = get_interval(opt)
% Get 95% conficence
% clc
% clear all
% close all



% figure(2222)
%     for i=1:size(Seq2,2)
%         subplot(3,ceil((size(Seq2,2)+1)/3),i)
%         if i <= size(Seq2,2)-2
%             plot(squeeze(Seq2(:,i,:)))
%         else
%             semilogy(squeeze(Seq2(:,i,:)))
%         end
%         set(gca,'xlim',[0 size(Seq2,1)*1.2])
%         set(gca,'fontsize',12)
%     end
%  
%     figure(1111)
%     for i=1:size(Seq1,2)
%         subplot(3,ceil((size(Seq1,2)+1)/3),i)
%         if i <= size(Seq1,2)-2
%             plot(squeeze(Seq1(:,i,:)))
%         else
%             semilogy(squeeze(Seq1(:,i,:)))
%         end
%         set(gca,'xlim',[0 size(Seq1,1)*1.2])
%         set(gca,'fontsize',12)
%     end


    load TEMP_RES.mat Sequences

    Seq = Sequences;


    index0=find(squeeze(Seq(:,end-1,1)),1,'first');
    index1=find(squeeze(Seq(:,1,1)),1,'last');

    Seq =  Seq(index0:index1,:,:);

    clear Sequences index0 index1


    tail = 150000*3;

    for i = 1:size(Seq,2)
            rSeq(:,i) = reshape([Seq(:,i,1)'; Seq(:,i,2)'; Seq(:,i,3)'], [], 1)';
    end


    rSeq = rSeq(end-tail:end,:);

    
figure()
for i = 1:11
    subplot(3,4,i)
    boxplot(rSeq(:,i))
end


%%
    bestset = rSeq(rSeq(:,end-1)==max(rSeq(:,end-1)), 1:end-2);
    bestset = bestset(end,:);

    if strcmp(opt,"best")==1
        x = bestset;
    elseif strcmp(opt,"rand")==1
        for i = 1:size(parmat1,2)
            k = randi([1 size(parmat,1)]);
            x(i) = parmat(k,i);
        end    
    elseif strcmp(opt,"95pct")==1
        for i = 1:size(rSeq,2)-2
            temp = sort(rSeq(:,i));
            pctmin(i)   = prctile(rSeq(:,i), 2.5);
            pctmax(i)   = prctile(rSeq(:,i), 97.5);
            parmat(:,i) = temp(round(tail*0.025):round(tail*0.975));
        end
        x(:,1) = pctmin; x(:,2) = pctmax;   
    elseif strcmp(opt,"50pct")==1
        for i = 1:size(rSeq,2)-2
            temp = sort(rSeq(:,i));
            pctmin(i)   = prctile(rSeq(:,i), 25);
            pctmax(i)   = prctile(rSeq(:,i), 75);
            parmat(:,i) = temp(round(tail*0.25):round(tail*0.75));
        end
        x(:,1) = pctmin; x(:,2) = pctmax;   
    end

end
