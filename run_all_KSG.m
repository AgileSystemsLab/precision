% KSG Wrapper:
% Run KSG method on all data in standard way
rng('shuffle')
nmoths = 7;
nmuscle = 10;
knn = 4;
noise = [0, logspace(log10(0.05), log10(6), 120)];
repeats = 150;
% Preallocate
MI = cell(nmoths, nmuscle);
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