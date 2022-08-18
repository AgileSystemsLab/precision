function [countofspikekx, countofspikenx, Knew, entropyvec] = conditionalspikesfunc(torquewordcolumn, words, j, g)

countofspikekx =[];
countofspikenx =[];
Knew = [];
%separating out words
wordswithmask = words;
entropyvec = [];
for x = min(torquewordcolumn):max(torquewordcolumn)
    mask1 = torquewordcolumn == x;
    wordswithmask = words(mask1);
    
    newprobdist = [];
      
    
    

    nx = [];%different counts in bins
    K = max(wordswithmask)-min(wordswithmask); %number of different counts
    kxword = []; %counts in different counts
    for u = min(wordswithmask):max(wordswithmask)
        kxword = [kxword, sum(wordswithmask == u)];
        newprobdist = [newprobdist, kxword(end)/length(wordswithmask)];% for direct calc

    end
    %direct method, not important for the estimator
    entropybyhand = 0; 
    for i = 1:length(newprobdist)
        if newprobdist(i) ~= 0
            entropybyhand = entropybyhand+newprobdist(i)*log(newprobdist(i));%/log(2);
        end
    end

    entropyvec = [entropyvec, entropybyhand];
    
    countofvectors = kxword;
    kxword = [];
    vectorofcounts = min(countofvectors):max(countofvectors);
    for u = vectorofcounts
        nx = [nx, sum(countofvectors == u)];
        kxword = [kxword, sum(countofvectors(countofvectors == u))];
        
    end
    mask = nx == 0;
    nx(mask) = [];
    kxword(mask) = [];
    mask = (kxword == 0);
    nx(mask) = [];
    kxword(mask) = [];
    %K = length(wordswithmask);
    K = sum(kxword);
    countofspikekx = [countofspikekx, {nx}];
    countofspikenx = [countofspikenx, {kxword}];
    Knew = [Knew, K];%we throw out this Knew, due to our changes in the way we're doing K
    
end
end

