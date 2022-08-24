function [S_nsbwordvec, dS_nsbwordvec, S_ml1wordvec, conditionalentropyvec, conditionaldS_nsbvec, conditionalS_ml1vec] = NSB_precision(X, Y, nspikingbins, ntorquebins)
    arguments
        X (:,:) double                  % X - array of spike times, each row is a wingbeat
        Y (:,:) double                  % Y - array of ouput variables, likely PC1,PC2 of yaw torque. Each row is a wb
        nspikingbins (1,1) double = 70  % Maximum number of word bins
        ntorquebins (1,1) double = 2;   % Maximum number of torque bins
    end

    % If NSB tools not on path, add to path 
    % (dumb and assumes folder within dir of this function)
    % (Checks entire path, could be way more efficient but ¯\_(ツ)_/¯ )
    path_cell = regexp(path, pathsep, 'split');
    function_dir = fileparts(mfilename('fullpath'));
    if ~any(strcmpi(fullfile(function_dir, 'nsb-entropy-code'), path_cell))
        addpath(fullfile(function_dir, 'nsb-entropy-code'))
    end
    tic

    
    %This next section sets up all the different data storage vectors
    conditionalentropyvec= zeros(1, nspikingbins);
    conditionalS_ml1vec = zeros(1, nspikingbins);
    conditionaldS_nsbvec = zeros(1, nspikingbins);
    
    errorcodevec = zeros(1, nspikingbins);
    S_nsbwordvec = zeros(1, nspikingbins);
    S_ml1wordvec = zeros(1, nspikingbins);
    dS_nsbwordvec = zeros(1, nspikingbins);
    
    binsizevec = [];

    count = 1;
    %Loops through different number of bins/bin sizes (g)
    for g = 1:nspikingbins
        
        numofbins = g;
        %determines the maximum time difference across the entire data set
        spiketimerange = max(max(X))-min(min(X));
        onebin = spiketimerange./numofbins;
        spiketimerange = [];
        binsizevec = [binsizevec, onebin];
        for x = 1:numofbins
            spiketimerange = [spiketimerange, x*onebin+min(min(X))];                
        end
        
        %creating array with each bin having the number of spikes, IOW creating an
        %array of words
        wordarray = [];
        
        for x = 1:length(X(:, 1))
            newwords = [];
            for y = 1:numofbins
                newwords = [newwords, sum(X(x,:) <= spiketimerange(y))];
            end
            wordarray = [wordarray; newwords];
        end
        newwords = [];
        holdup = wordarray(:, 1);
        x = length(wordarray(1, :));
        while x > 1
            newwords = [ wordarray(:, x)-wordarray(:, x-1), newwords];
            x = x-1;
        end
        
        
        wordarray = [holdup, newwords];
        %Find the number of different words in the set
        [newwords, ia, ~] = (unique(wordarray, 'rows'));
        newwords = flipud(newwords);
        ia = flipud(ia);
        wordname = [1:length(ia)];
        words = zeros(length(wordarray(:,1)), 1);
        %turns the multidimensional array of words into a #wingstrokesxword array
        %where each array is a wingstroke with the designated word
        for o = 1:length(wordname)
            mask = ismember(wordarray,newwords(o,:),'rows');
            words(mask) = o;
        end
        
        nx = [];%different counts in bins
        kxword = []; %counts in different counts
        for u = min(words):max(words)
            kxword = [kxword, sum(words == u)];
        end
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
        
        K = 24360;%large K, feel free to vary, I did not observe any changes

        
        [S_nsbword, dS_nsbword, ~, ~, ~, S_ml1,errcodeword] = find_nsb_entropy (nx, kxword, K, .1, 1); %input names are reversed here
        errorcodevec(count) = errcodeword;
        S_nsbwordvec(count) = S_nsbword;
        S_ml1wordvec(count) = S_ml1;
        dS_nsbwordvec(count) = dS_nsbword;
        
        for j = ntorquebins %this is where the word spikes are split into torque words
            tempstdvec =[];
            [probdist, torquewordcolumn] = torquebreakups(j, Y);%this function just splits up the torque into j bins and then uses that to create torque words
            [countofspikekx, countofspikenx, ~, ~] = conditionalspikesfunc(torquewordcolumn, words, j, g); % this function identifies the spiking words that identify with the torquewords
            
            conditionalentropy = 0;
            conditionalentropysml1 = 0;
            for f = 1:length(countofspikekx)% this is where the estimation happens
                 if isequal(countofspikekx{f} == 1, ones(length(countofspikekx{f}))) && isequal(countofspikenx{f} == 1, ones(length(countofspikenx{f})))
                     S_nsbword = 0;
                     S_ml1 = 0;
                 else
                     Knew = 24360;
                    [S_nsbword, dS_nsbword, ~, ~, ~, S_ml1,errcodeword] = find_nsb_entropy (countofspikekx{f}, countofspikenx{f}, Knew, .1, 1);
                 end
                if errcodeword ~= 0
                   posd = 8; 
                end
                tempstdvec = [tempstdvec, ((dS_nsbword)^2)*probdist(f)];%standard deviation calculation
                conditionalentropy = conditionalentropy+S_nsbword*probdist(f);
                conditionalentropysml1 = conditionalentropysml1+S_ml1*probdist(f);
                
            end
            conditionalentropyvec(count) = conditionalentropy;%storing everything
            conditionalS_ml1vec(count) = conditionalentropysml1;
            conditionaldS_nsbvec(count) = sqrt(sum(tempstdvec));

        end
        count = count+1;
    end
    
    %converting everything to bits
    S_nsbwordvec = S_nsbwordvec./log(2);
    dS_nsbwordvec = dS_nsbwordvec./log(2);
    S_ml1wordvec = S_ml1wordvec./log(2);
    conditionalentropyvec = conditionalentropyvec./log(2);
    conditionaldS_nsbvec = conditionaldS_nsbvec./log(2);
    conditionalS_ml1vec = conditionalS_ml1vec./log(2);
    toc
end