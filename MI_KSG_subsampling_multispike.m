function [MIs] = MI_KSG_subsampling_multispike(X, Y, knn, split_sizes, do_plot)
    arguments 
        X (:,:) double
        Y (:,:) double
        knn (1,1) double = 4
        split_sizes (1,:) double = [1,2,3,4,5]
        do_plot (1,1) logical = false
    end
    % Preallocate, setup
    n = length(split_sizes);
    MIs = cell(n,2);
    if do_plot
        means = zeros(n,1);
        errorbars = zeros(n,1);
    end
    % If KSG tools not on path, add to path 
    % (dumb and assumes folder within dir of this function)
    % (Checks entire path, could be way more efficient but ¯\_(ツ)_/¯ )
    path_cell = regexp(path, pathsep, 'split');
    function_dir = fileparts(mfilename('fullpath'));
    if ~any(strcmpi(fullfile(function_dir, 'ContinuousMIEstimation'), path_cell))
        addpath(fullfile(function_dir, 'ContinuousMIEstimation'))
    end
    % Loop over data splits
    for i = 1:n
        % Permute data, assign at random to one of the subsets
        a = randperm(length(X));
        nsplits = split_sizes(i);
        l = round(linspace(0, length(X), nsplits+1));
        % preallocate, setup MI cell array
        MIs{i,1} = i;
        MIs{i,2} = zeros(nsplits, 1);
        % Loop over each subset in this split
        for j = 1:nsplits
            % X and Y data to use in this subset
            Xsub = X(a(l(j) + 1 : l(j+1)), :);
            Ysub = Y(a(l(j) + 1 : l(j+1)), :);
            % Vector of how many spikes per wingbeat from X
            Nspike = sum(~isnan(Xsub), 2);
            % Calculate probability of each possible number of spikes
            unq = unique(Nspike);
            probs = arrayfun(@(x) sum(Nspike==x) / length(Nspike), unq);
            probsmap = containers.Map(unq, probs);
            % Loop over number of spikes in wb
            % (Skip 0 spike group, because 0 spikes means no information)
            for k = unq(unq~=0)'
                % Get indices of data with this many spikes
                inds = Nspike == k;
                % Call mutual information estimator
                if sum(inds) > k && knn < sum(inds)
                    mi = MIxnyn_matlab(Xsub(inds, 1:k), Ysub(inds,:), knn);
                else
                    mi = 0;
                end
                % Convert from nats to bits
                mi = mi / log(2);
                % Adjust based on probability of this many spikes
                MIs{i,2}(j) = MIs{i,2}(j) + mi * probsmap(k);
            end
        end
        % Get mean and std of subsamples if making plot
        if do_plot
            means(i) = mean(MIs{i,2});
            errorbars(i) = std(MIs{i,2});
        end
    end
    if do_plot
        figure
        errorbar(split_sizes ./ length(X), means, errorbars)
        xlabel('1/N')
        ylabel('Mutual Information, bits')
        title(['Dependence of Mutual Information on the data set size for k = ', num2str(knn)])
    end
end