function [conditionalentropyvec_bias, conditionalS_ml1vec_bias, conditionaldS_nsbvec_bias, conditionvariance_bias] = NSB_precision_bias(X, Y, nspikingbins, ntorquebins, repeats)
    arguments
        X (:,:) double                  % X - array of spike times, each row is a wingbeat
        Y (:,:) double                  % Y - array of ouput variables, likely PC1,PC2 of yaw torque. Each row is a wb
        nspikingbins (1,1) double = 70  % Maximum number of word bins
        ntorquebins (1,1) double = 2;   % Maximum number of torque bins
        repeats (1,1) double = 5;       % How many times to repeat to estimate bias
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
    conditionalentropyvec_bias = zeros(nspikingbins, repeats);
    conditionalS_ml1vec_bias = zeros(nspikingbins, repeats);
    conditionaldS_nsbvec_bias = zeros(nspikingbins, repeats);

    
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
        wordname = 1:length(ia);
        words = zeros(length(wordarray(:,1)), 1);
        %turns the multidimensional array of words into a #wingstrokesxword array
        %where each array is a wingstroke with the designated word
        for o = 1:length(wordname)
            mask = ismember(wordarray,newwords(o,:),'rows');
            words(mask) = o;
        end
        
        for j = ntorquebins %this is where the word spikes are split into torque words
            for jj = 1:repeats
                tempstdvec = 0;
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
                        [S_nsbword, dS_nsbword, ~, ~, ~, S_ml1, ~] = find_nsb_entropy(countofspikekx{f}, countofspikenx{f}, Knew, .1, 1);
                     end
                    tempstdvec = tempstdvec + ((dS_nsbword)^2)*probdist(f); %standard deviation calculation
                    conditionalentropy = conditionalentropy+S_nsbword*probdist(f);
                    conditionalentropysml1 = conditionalentropysml1+S_ml1*probdist(f);
                    
                end
                % Store everything
                conditionalentropyvec_bias(count, jj) = conditionalentropy;
                conditionalS_ml1vec_bias(count, jj) = conditionalentropysml1;
                conditionaldS_nsbvec_bias(count, jj) = sqrt(tempstdvec);
            end
        end
        count = count+1;
    end
    
    %converting everything to bits
    conditionalentropyvec_bias = conditionalentropyvec_bias./log(2);
    conditionalS_ml1vec_bias = conditionalS_ml1vec_bias./log(2);
    conditionaldS_nsbvec_bias = conditionaldS_nsbvec_bias./log(2);
    conditionvariance_bias = conditionaldS_nsbvec_bias.^2;
    toc
end