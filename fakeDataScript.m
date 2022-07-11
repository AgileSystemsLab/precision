rng('shuffle')
% Moth and muscle to focus on
moth = '1';
muscle = 'RAX';
noise = (0:.05:6);
repeats = 150;
n = 10; % How many different levels of rescaling to use

%---- Load data
load(fullfile('Data',['Moth',num2str(moth),'_MIdata.mat']))
% Rename so it's easier to write, get some useful quantities out
X = time_data.([muscle,'strokes']);
Y = Tz_WSd;
Nspike = sum(~isnan(X), 2);

% %---- Fake data: Totally uncorrelated
% % Same resolution, range, and # of spikes as real data, but uniformly distributed
% fake_uncorr = nan(size(X));
% % Shift so min resolution is 1, allowing rand generation using ints
% minmax = round([min(X, [], 'all'), max(X, [], 'all')] * 10); 
% for i = 1:size(X,1)
%     fake_uncorr(i,1:Nspike(i)) = randi(minmax, 1, Nspike(i));
% end
% fake_uncorr = fake_uncorr / 10;

%---- Direct connected fake data compared to real
% Setup parameters and preallocate
data_fraction = linspace(0.4, 1, n);
MI_direct = cell(1,n);
MI_direct(:) = {zeros(length(noise), repeats)};
fake_X = cell(1,n);
% fake_X(:) = {nan(size(X))};
% Loop over rescaling, create fake data, run precision estimation
for i = 1:n
    len_subsample = round(data_fraction(i) * length(Y));
    subsample = randperm(length(Y), len_subsample);
    fake_X{i} = X(subsample,:);
    MI_direct{i} = KSG_precision(fake_X{i}, Y(subsample,:), 4, 20, noise);
end


%---- Run precision estimation on real and fake data and plot each
MI_real = KSG_precision(X, Y, 4, repeats, noise, true);
title('Real data')
% MI_uncorr = KSG_precision(fake_uncorr, Y, 4, repeats, noise, true);
% title('Uncorrelated fake data')

%---- Deterministically linked fake data at different levels of noise corruption
figure()
hold on
cols = copper(n);
for i = 1:n
    signal_amp = max(fake_X{i}, [], 'all') - min(fake_X{i}, [], 'all');
    rescale = 1 / mean(MI_direct{i}(1,:));
%     rescale = 1;
    plot(log10(noise), rescale * mean(MI_direct{i}, 2), 'color', cols(i,:));
end
title('motor PCs directly to spike times with varying m')
xlabel('log10( noise added / signal range )')
ylabel('normalized information (1 at zero noise added)')
colormap(copper)
cbh = colorbar;
set(cbh,'YTick',linspace(0,1,n))
set(cbh,'YTickLabel', num2str(data_fraction.'))

%% Distribution of kth nearest neighbor distances as noise is added
ncheck = 5;
figure
hold on
xlims = zeros(ncheck, 2);
ax = gobjects(ncheck, 1);
noiseinds = round(linspace(1,ncheck,ncheck)/ncheck*length(noise));
for ii = 1:ncheck
    noiseamp = noise(noiseinds(ii));
    noise_X = X(:,1) + noiseamp * rand(length(X), 1);
    kdist = zeros(length(X),1);
    for i = 1:length(X)
        dist = mink(noise_X(i,1) - noise_X(:,1), 4);
        kdist(i) = abs(dist(1));
    end
    ax(ii) = subaxis(ncheck, 1 , ii, 'SpacingVert', 0);
    histogram(kdist, 40, 'Normalization', 'pdf', 'EdgeColor', 'none')
    xline(median(kdist, 'omitnan'), 'k')
    text(0.6, 0.5, ['r = ', num2str(noiseamp)], 'units', 'normalized')
    xlims(ii,:) = get(gca, 'Xlim');
end
% Link axes so X axis is shared
linkaxes(ax)
xlabel(ax(end), 'abs() of k=4th nearest neighbor distances')
ylabel(ax(round(length(ax)/2)), 'Probability')


