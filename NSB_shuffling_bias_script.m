% NSB Shuffling Bias Script
addpath('nsb-entropy-code')
numoftorquebins = 5;%maximum number of torque bins
numofspikingbins = 70;%maximum number of word bins
repeats = 3;
numofmoths = 7;

%This next section sets up all the different data storage vectors
conditionalentropyvec = zeros(numofmoths*10, numofspikingbins, numoftorquebins, repeats);
conditionalentropyerrorcodevec = zeros(numofmoths*10, numofspikingbins, numoftorquebins, repeats);
conditionalS_ml1vec = zeros(numofmoths*10, numofspikingbins, numoftorquebins, repeats);
conditionaldS_nsbvec = zeros(numofmoths*10, numofspikingbins, numoftorquebins, repeats);
rangevec = zeros(numofmoths*10, numofspikingbins);

othercount = 1;
binsizevec = [];
for i = 1:numofmoths %same data extraction process
    load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
    mothname = ['Moth', num2str(i)];
    fields = fieldnames(time_data);
    for fn=fields'
        fn;
        %fn = {'LDLMstrokes'} %Change this to a specific muscle to debug
        %that muscle
        %# since +fn is a 1-by-1 cell array, you still need to index into it, unfortunately
        timedata_data = time_data.(fn{1});
        count = 1;
        musclename = fn{1};
        musclename = musclename(1:end-7);
        
        arrayofwords = [];
        binvec = [];
        
        conditionalentropy = 0;
        %Loops through different number of bins/bin sizes (g)
        for g = 1:numofspikingbins
            
            numofbins = g;
            %determines the maximum time difference across the entire data set
            spiketimerange = max(max(timedata_data))-min(min(timedata_data));
            onebin = spiketimerange./numofbins;
            spiketimerange = [];
            binsizevec = [binsizevec, onebin];
            rangevec(othercount, count) = onebin;
            %designating the range of forces
            
            for x = 1:numofbins
                spiketimerange = [spiketimerange, x*onebin+min(min(timedata_data))];
            end
            
            numofelements = sum(sum(~isnan(timedata_data)));
            numinbin = [];
            
            %creating array with each bin having the number of spikes, IOW creating an
            %array of words
            wordarray = [];
            
            for x = 1:length(timedata_data(:, 1))
                newwords = [];
                for y = 1:numofbins
                    newwords = [newwords, sum(timedata_data(x,:) <= spiketimerange(y))];
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
            [newwords, ia, ic] = (unique(wordarray, 'rows'));
            newwords = flipud(newwords);
            ia = flipud(ia);
            numwords = flipud(accumarray(ic, 1));
            
            wordpmf = numwords./sum(numwords);
            
            words = zeros(length(wordarray(:,1)), 1);
            
            %turns the multidimensional array of words into a #wingstrokesxword array
            %where each array is a wingstroke with the designated word
            wordnames = rand(1, length(ia));
            for o = 1:length(ia)
                
                mask = ismember(wordarray,newwords(o,:),'rows');
                
                words(mask) = o;
                
            end
            
            nx = [];%different counts in bins
%             K = max(words)-min(words); %number of different counts
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

%             [S_nsbword, dS_nsbword, S_clword, dS_clword, xi_clword, S_ml1,errcodeword] = find_nsb_entropy (nx, kxword, K, .1, 1); %input names are reversed here
%             errorcodevec(othercount, count) = errcodeword;
%             S_nsbwordvec(othercount, count) = S_nsbword;
%             S_ml1wordvec(othercount, count) = S_ml1;
%             dS_nsbwordvec(othercount, count) = dS_nsbword;
            
            for j = 2:numoftorquebins %this is where the word spikes are split into torque words
                for jj = 1:repeats
                    tempstdvec = [];
                    [probdist, torquewordcolumn] = torquebreakups(j, Tz_WSd(randperm(length(Tz_WSd)),:));%this function just splits up the torque into j bins and then uses that to create torque words
    %                 [probdist, torquewordcolumn] = torquebreakups1columnonlyplease(j, Tz_WSd(:, 1));%this does only first pc column
                    [countofspikekx, countofspikenx, Knew, ~] = conditionalspikesfunc(torquewordcolumn, words, j, g); % this function identifies the spiking words that identify with the torquewords
                    
                    conditionalentropy = 0;
                    conditionalentropysml1 = 0;
                    for f = 1:length(countofspikekx)% this is where the estimation happens
                         if isequal(countofspikekx{f} == 1, ones(length(countofspikekx{f}))) && isequal(countofspikenx{f} == 1, ones(length(countofspikenx{f})))
                             S_nsbword = 0;
                             S_ml1 = 0;
                         else
                             Knew = 24360;
                            [S_nsbword, dS_nsbword, S_clword, dS_clword, xi_clword, S_ml1, ~] = find_nsb_entropy(countofspikekx{f}, countofspikenx{f}, Knew, .1, 1);
                         end
                        tempstdvec = [tempstdvec, ((dS_nsbword)^2)*probdist(f)];%standard deviation calculation
                        conditionalentropy = conditionalentropy+S_nsbword*probdist(f);
                        conditionalentropysml1 = conditionalentropysml1+S_ml1*probdist(f);
                        
                    end
                    conditionalentropyvec(othercount, count, j, jj) = conditionalentropy;%storing everything
                    conditionalS_ml1vec(othercount, count, j, jj) = conditionalentropysml1;
                    conditionaldS_nsbvec(othercount, count, j, jj) = sqrt(sum(tempstdvec));
                end
            end
            count = count+1;
        end
        othercount = othercount +1;
    end
end

%converting everything to bits
conditionalentropyvec_bias = conditionalentropyvec./log(2);
conditionalS_ml1vec_bias = conditionalS_ml1vec./log(2);
conditionaldS_nsbvec_bias = conditionaldS_nsbvec./log(2);
conditionvariance_bias = conditionaldS_nsbvec.^2;
% Save results
save('NSB_bias_data.mat', ...
    'conditionalentropyvec_bias', 'conditionalS_ml1vec_bias', ...
    'conditionaldS_nsbvec_bias', 'conditionvariance_bias');
