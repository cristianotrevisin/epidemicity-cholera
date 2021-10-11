function [x1,x2] = get95(opt,rSeq1, rSeq2)
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
tail2 = 25000*3;

    tail1 = 1000000*3;



bestset1 = rSeq1(rSeq1(:,end-1)==max(rSeq1(:,end-1)), 1:end-2);
bestset1 = bestset1(end,:);

bestset2 = rSeq2(rSeq2(:,end-1)==max(rSeq2(:,end-1)), 1:end-2);
bestset2 = bestset2(end,:);

for i = 1:size(rSeq2,2)-2
    temp = sort(rSeq2(:,i));
    pctmin2(i)   = prctile(rSeq2(:,i), 5);
    pctmax2(i)   = prctile(rSeq2(:,i), 95);
    parmat2(:,i) = temp(round(tail2*0.05):round(tail2*0.95));
end

for i = 1:size(rSeq1,2)-2
    temp = sort(rSeq1(:,i));
    pctmin1(i)   = prctile(rSeq1(:,i), 5);
    pctmax1(i)   = prctile(rSeq1(:,i), 95);
    parmat1(:,i) = temp(round(tail1*0.05):round(tail1*0.95));
end



if strcmp(opt,"best")==1
    x1 = bestset1;
    x2 = bestset2;
elseif strcmp(opt,"rand")==1
    for i = 1:size(parmat1,2)
        k1 = randi([1 size(parmat1,1)]);
        x1(i) = parmat1(k1,i);
    end
    for i = 1:size(parmat2,2)
        k2 = randi([1 size(parmat2,1)]);
        x2(i) = parmat2(k2,i);
    end
    
elseif strcmp(opt,"95pct")==1
    x1(:,1) = pctmin1; x1(:,2) = pctmax1;
    x2(:,1) = pctmin2; x2(:,2) = pctmax2;   
end

end
