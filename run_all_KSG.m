% KSG Wrapper. Run KSG method on all data
rng('shuffle')
nmoths = 7;
nmuscles = 10;
noise = [0, logspace(log10(0.05), log10(6), 120)];
repeats = 150;

run_standard = false;
run_varyk = false;
run_getcountinfo = true;

%% Run KSG method on all data in standard way
if run_standard
    knn = 4;
    % Preallocate
    MI = cell(nmoths, nmuscles);
    MI(:) = {zeros(length(noise), repeats)};
    % Loop over moths
    for i = 1:nmoths
        % Load data
        load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
        fields = fieldnames(time_data);
        disp(['Moth ',num2str(i)])
        % Loop over muscles
        for j = 1:length(fields)
            disp(fields{j}(1:end-7))
            % Run KSG at noise levels to find precision
            MI{i,j} = KSG_precision(time_data.(fields{j}), Tz_WSd, knn, repeats, noise);
        end
    end
    % Save results
    save('KSG_data.mat', 'MI', 'noise', 'repeats', 'knn')
end

%% Get spike count information
if run_getcountinfo
    knn = 4;
    MI_count = zeros(nmoths, nmuscles);
    % If KSG tools not on path, add to path 
    % (dumb and assumes folder within dir of this function)
    % (Checks entire path, could be way more efficient but ¯\_(ツ)_/¯ )
    path_cell = regexp(path, pathsep, 'split');
    function_dir = fileparts(mfilename('fullpath'));
    if ~any(strcmpi(fullfile(function_dir, 'ContinuousMIEstimation'), path_cell))
        addpath(fullfile(function_dir, 'ContinuousMIEstimation'))
    end
    % Loop over moths
    for i = 1:nmoths
        % Load data
        load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
        fields = fieldnames(time_data);
        disp(['Moth ',num2str(i)])
        % Loop over muscles
        for j = 1:length(fields)
            disp(fields{j}(1:end-7))
            % Get count information 
            MI_count(i,j) = MIxnyn_matlab(sum(~isnan(time_data.(fields{j})),2), Tz_WSd, knn);
        end
    end
    % Save results
    save('KSG_data_count.mat', 'MI_count', 'knn')
end



%% Run KSG method on all data while varying k
if run_varyk
    knn_vec = [2,3,4,5,6,7];
    % Preallocate
    MI_varyk = cell(nmoths, nmuscles, length(knn_vec));
    % Loop over levels of k
    for k = 1:length(knn_vec)
        % Loop over moths
        for i = 1:nmoths
            % Load data
            load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
            fields = fieldnames(time_data);
            disp(['Moth ',num2str(i)])
            % Loop over muscles
            for j = 1:length(fields)
                disp(fields{j}(1:end-7))
                % Run KSG at noise levels to find precision
                MI_varyk{i,j,k} = KSG_precision(time_data.(fields{j}), Tz_WSd, knn_vec(k), repeats, noise);
            end
        end
    end
    % Save results
    save('KSG_data_varyk.mat', 'MI_varyk', 'noise', 'repeats', 'knn_vec')
end