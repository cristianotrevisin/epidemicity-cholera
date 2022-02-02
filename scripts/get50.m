function [x1,x2] = get95(opt,rSeq1, rSeq2)
% GET95 Gets 95% conficence interval
% Takes as arguments:
%       - opt: the option (best for best set; 95pct for 50% CI; rand for
%       random generation)
%       - rSeq1: processed sequence from the first MC
%       - rSeq2: processed sequence from the second MC
%
% Returns the following outputs:
%       - x1, x2 the requested generated elements.

    tail2 = 20000*3;
    tail1 = 1000000*3;

    bestset1 = rSeq1(rSeq1(:,end-1)==max(rSeq1(:,end-1)), 1:end-2);
    bestset1 = bestset1(end,:);

    bestset2 = rSeq2(rSeq2(:,end-1)==max(rSeq2(:,end-1)), 1:end-2);
    bestset2 = bestset2(end,:);

    for i = 1:size(rSeq2,2)-2
        temp = sort(rSeq2(:,i));
        pctmin2(i)   = prctile(rSeq2(:,i), 25);
        pctmax2(i)   = prctile(rSeq2(:,i), 75);
        parmat2(:,i) = temp(round(tail2*0.25):round(tail2*0.75));
    end

    for i = 1:size(rSeq1,2)-2
        temp = sort(rSeq1(:,i));
        pctmin1(i)   = prctile(rSeq1(:,i), 25);
        pctmax1(i)   = prctile(rSeq1(:,i), 75);
        parmat1(:,i) = temp(round(tail1*0.25):round(tail1*0.75));
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
