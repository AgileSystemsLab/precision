rng('shuffle')
% Moth and muscle to focus on
moth = '1';
muscle = 'RAX';
knn = 4;
noise = (0:.05:6);
repeats = 150;
n = 3; % How many different levels of rescaling to use

%---- Load data
load(fullfile('Data',['Moth',num2str(moth),'_MIdata.mat']))
% Rename so it's easier to write, get some useful quantities out
X = time_data.([muscle,'strokes']);
Y = Tz_WSd;
Nspike = sum(~isnan(X), 2);
unq = unique(Nspike);

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
    MI_direct{i} = KSG_precision(fake_X{i}, Y(subsample,:), knn, 100, noise);
end


%---- Run precision estimation on real and fake data and plot each
MI_real = KSG_precision(X, Y, knn, repeats, noise, true);
title('Real data')
% MI_uncorr = KSG_precision(fake_uncorr, Y, knn, repeats, noise, true);
% title('Uncorrelated fake data')

%---- Deterministically linked fake data at different sample sizes
figure()
hold on
cols = copper(n);
for i = 1:n
    signal_amp = max(fake_X{i}, [], 'all') - min(fake_X{i}, [], 'all');
    rescale = 1 / mean(MI_direct{i}(1,:));
%     rescale = 1;
    plot(log10(noise), rescale * mean(MI_direct{i}, 2), 'color', cols(i,:));
end
title('motor PCs directly to spike times with varying sample size')
xlabel('log10( noise added / signal range )')
ylabel('normalized information (1 at zero noise added)')
colormap(copper)
cbh = colorbar;
set(cbh,'YTick',linspace(0,1,n))
set(cbh,'YTickLabel', num2str(data_fraction.'))



%% Distribution of kth nearest neighbor distances as noise is added
% Levels of noise to check
ncheck = 5;
kdist = cell(ncheck, length(unq(unq~=0)));
% Loop over number of spikes in a wingbeat
for jj = unq(unq~=0)'
    useX = X(Nspike==jj, 1:jj);
    useY = Y(Nspike==jj, :);
    % Continue only if enough wingbeats with this many spikes
    if size(useX, 1) >= knn
        % Preallocate, setup figure
        figure
        hold on
        xlims = zeros(ncheck, 2);
        ax = gobjects(ncheck, 1);
        noiseinds = linspace(1, length(noise), ncheck);
        % Loop over noise levels
        for ii = 1:ncheck
            noiseamp = noise(noiseinds(ii));
            noiseX = useX + noiseamp * rand(size(useX));
            kdist{ii,jj} = nan(size(useX,1),1);
            % Loop over each point in X, get knn distances
            for i = 1:size(useX,1)
                % Get indices of kth nearest neighbors, correct indices (X and Y)
                [idx, distx] = knnsearch(noiseX([1:(i-1), (i+1):end], :), noiseX(i,:), 'K', knn);
                idx(idx>i) = idx(idx>i) + 1;
                [idy, disty] = knnsearch(useY([1:(i-1), (i+1):end], :), useY(i,:), 'K', knn);
                idy(idy>i) = idy(idy>i) + 1;
                % Select max distance
                kdist{ii,jj}(i) = max(distx(knn), disty(knn));
            end
            % Make plot
            ax(ii) = subaxis(ncheck, 1 , ii, 'SpacingVert', 0);
            histogram(log10(kdist{ii,jj}), 40, 'Normalization', 'pdf', 'EdgeColor', 'none')
            xline(median(log10(kdist{ii,jj}), 'omitnan'), 'k')
            text(0.6, 0.5, ['r = ', num2str(noiseamp)], 'units', 'normalized')
            xlims(ii,:) = get(gca, 'Xlim');
        end
        % Link axes so X axis is shared, clean up X ticks
        xtickspace = linspace(10^min(xlims(:,1)), 10^max(xlims(:,2)), 8);
        for ii = 1:ncheck
            ax(ii).set('xlim', [min(xlims(:,1)), max(xlims(:,2))]);
            ax(ii).set('xtick', log10(xtickspace));
            ax(ii).set('xticklabels', []);
        end
        ax(end).set('xticklabels', num2cell(round(10.^get(ax(end),'XTick'), 2)));
        xlabel(ax(end), ...
            ['k=',num2str(knn),'$^{th}$ N.N. distance in space $\|z-z`\| = max\{\|x-x`\|, \|y-y`\|\}$'], ...
            "interpreter", "latex")
        ylabel(ax(round(length(ax)/2)), 'Probability')
        title(ax(1), [num2str(jj), ' spikes in wingbeat'])
    end
end


%% How MI and precision improve with transformations
% From Kraskov, replacing data with rank (index of order in total dataset)
% can improve estimations, as MI should be invariant but transform gives
% uniform density
% Also from Kraskov, transforms can be a good idea in cases where
% distributions are skewed, non-symmetric, or rough. log transform should
% lead to symmetric, normal knn distances (as seen above, where distance
% distribution is lognormal)
[~,~,I] = unique(X + 1e-8 * rand(size(X)));

xrange = range(X, 'all');
logxrange = range(log(X - min(X) + 0.1), 'all');

% ranknoise = noise / range(X) * range(I);
lognoise = noise / xrange * logxrange;
% MI_rank = KSG_precision(I, Y, knn, repeats, ranknoise, true);
MI_log = KSG_precision(log(X - min(X) + 0.1), log(Y - min(Y) + 0.1), knn, repeats, lognoise);

figure
mseb(log10(noise/xrange), mean(MI_real,2), std(MI_real,0,2)', struct('col',{{'b'}}))
% mseb(log10(ranknoise/range(I)), mean(MI_rank,2), std(MI_rank,0,2)', struct('col',{{'r'}}))
mseb(log10(lognoise/logxrange), mean(MI_log,2), std(MI_log,0,2)', struct('col',{{'g'}}))
xlabel('log10(noise amplitude / signal amplitude)')
ylabel('Mutual Information (bits)')

mean_MI_real = mean(MI_real, 2);
% mean_MI_rank = mean(MI_rank, 2);
mean_MI_log = mean(MI_log, 2);
real_p = noise(find(mean_MI_real < (mean_MI_real(1) - std(MI_real(1,:))), 1))
% rank_p = ranknoise(find(mean_MI_rank < (mean_MI_rank(1) - std(MI_rank(1,:))), 1)) / range(I) * range(X)
log_p = lognoise(find(mean_MI_log < (mean_MI_log(1) - std(MI_log(1,:))), 1)) / logxrange * xrange
