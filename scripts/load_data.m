function [rSeq1, rSeq2] = load_data();
% LOAD_DATA Loads the data from the Markov Chains.
%   Inputs: none.
%   Outputs:
%       - rSeq1: processed sequence from first Markov Chain;
%       - rSeq2: processed sequence from second Markov Chain;

    load ../markov_chains/TEMP_RES_2.mat Sequences
    Seq2 = Sequences;

    index0=find(squeeze(Seq2(:,end-1,1)),1,'first');
    index1=find(squeeze(Seq2(:,1,1)),1,'last');

    Seq2 =  Seq2(index0:index1,:,:);

    clear Sequences index0 index1


    load ../markov_chains/TEMP_RES_1.mat Sequences
    Seq1 = Sequences;


    index0=find(squeeze(Seq1(:,end-1,1)),1,'first');
    index1=find(squeeze(Seq1(:,1,1)),1,'last');

    Seq1 =  Seq1(index0:index1,:,:);

    clear Sequences index0 index1


    tail2 = 20000*3;

    tail1 = 1000000*3;

    for i = 1:size(Seq1,2)
        rSeq1(:,i) = reshape([Seq1(:,i,1)'; Seq1(:,i,2)'; Seq1(:,i,3)'], [], 1)';
    end

    for i = 1:size(Seq2,2)
        rSeq2(:,i) = reshape([Seq2(:,i,1)'; Seq2(:,i,2)'; Seq2(:,i,3)'], [], 1)';
    end

    rSeq1 = rSeq1(end-tail1:end,:);
    rSeq2 = rSeq2(end-tail2:end,:);
    
end