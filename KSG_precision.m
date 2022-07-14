% KSG Mutual Information Estimation of Precision via Uniform Noise Corruption
function [MI] = KSG_precision(X, Y, knn, repeats, noise, doplot, runparallel, dolog)
    arguments
        X (:,:) double                  % X - array of spike times, each row is a wingbeat
        Y (:,:) double                  % Y - array of ouput variables, likely PC1,PC2 of yaw torque. Each row is a wb
        knn (1,1) double = 4            % knn - number of nearest neighbors
        repeats (1,1) double = 150      % repeats - How many times to repeat MI estimator
        noise (1,121) double = (0:.05:6)% noise - uniform noise levels to apply to spike times
        doplot (1,1) logical = false    % whether or not to make a plot
        runparallel (1,1) logical = true% Whether or not to parallelize running over noise levels
        dolog (1,1) logical = false
    end
    tic
    % If KSG tools not on path, add to path 
    % (dumb and assumes folder within dir of this function)
    % (Checks entire path, could be way more efficient but ¯\_(ツ)_/¯ )
    path_cell = regexp(path, pathsep, 'split');
    function_dir = fileparts(mfilename('fullpath'));
    if ~any(strcmpi(fullfile(function_dir, 'ContinuousMIEstimation'), path_cell))
        addpath(fullfile(function_dir, 'ContinuousMIEstimation'))
    end
    % Create vector of how many spikes per wingbeat from X
    Nspike = sum(~isnan(X), 2);
    % Calculate probability of each possible number of spikes
    unq = unique(Nspike);
    probs = arrayfun(@(x) sum(Nspike==x)/length(Nspike), unq);
    probsmap = containers.Map(unq, probs);
    % Preallocate
    MI = zeros(length(noise), repeats);
    if dolog
        if isempty(gcp('nocreate'))
            parpool();
        end
        % Loop over noise levels
        parfor i = 1:length(noise)
            % Loop over repeats
            for j = 1:repeats
                % Loop over number of spikes in wb
                % Skip 0 spike group, because 0 spikes means no information
                for k = unq(unq~=0)'
                    % Get indices of data with this many spikes
                    inds = Nspike == k;
                    % Call mutual information estimator
                    if sum(inds) > k && knn < sum(inds)
                        mi = MIxnyn_matlab(...
                            log(X(inds, 1:k) + noise(i)*rand(sum(inds), k) - min(X(inds,1:k)) + 0.1), ...
                            log(Y(inds,:) - min(Y(inds,:)) + 0.1), knn, pwd);
                    else
                        mi = 0;
                    end
                    % Adjust based on probability of this many spikes
                    MI(i,j) = MI(i,j) + mi * probsmap(k);
                end
            end
        end
    else
        if runparallel
            if isempty(gcp('nocreate'))
                parpool();
            end
            % Loop over noise levels
            parfor i = 1:length(noise)
                % Loop over repeats
                for j = 1:repeats
                    % Loop over number of spikes in wb
                    % Skip 0 spike group, because 0 spikes means no information
                    for k = unq(unq~=0)'
                        % Get indices of data with this many spikes
                        inds = Nspike == k;
                        % Call mutual information estimator
                        if sum(inds) > k && knn < sum(inds)
                            mi = MIxnyn_matlab(X(inds, 1:k) + noise(i)*rand(sum(inds), k), Y(inds,:), knn, pwd);
                        else
                            mi = 0;
                        end
                        % Adjust based on probability of this many spikes
                        MI(i,j) = MI(i,j) + mi * probsmap(k);
                    end
                end
            end
        % Serial version
        else
            % Loop over noise levels
            for i = 1:length(noise)
                % Loop over repeats
                for j = 1:repeats
                    % Loop over number of spikes in wb
                    % Skip 0 spike group, because 0 spikes means no information
                    for k = unq(unq~=0)'
                        % Get indices of data with this many spikes
                        inds = Nspike == k;
                        % Call mutual information estimator
                        if sum(inds) > k && knn < sum(inds)
                            mi = MIxnyn_matlab(X(inds, 1:k) + noise(i)*rand(sum(inds), k), Y(inds,:), knn, pwd);
                        else
                            mi = 0;
                        end
                        % Adjust based on probability of this many spikes
                        MI(i,j) = MI(i,j) + mi * probsmap(k);
                    end
                end
            end
        end
    end
    % Convert from nats to bits
    MI = MI / log(2);
    % Make plot, if requested
    if doplot
        figure();
        mseb(log10(noise), mean(MI,2), std(MI, 0, 2)');
    end
    toc
end