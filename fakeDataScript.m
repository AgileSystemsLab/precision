rng('shuffle')
% Moth and muscle to focus on
moth = '1';
muscle = 'LBA';
knn = 4;
% noise = (0:.05:6);
noise = [0, logspace(log10(0.05), log10(6), 120)];
repeats = 150;
n = 15; 

%---- Load data
load(fullfile('Data',['Moth',num2str(moth),'_MIdata.mat']))
% Rename so it's easier to write, get some useful quantities out
X = time_data.([muscle,'strokes']);
Y = Tz_WSd;
Nspike = sum(~isnan(X), 2);
unq = unique(Nspike);

% %---- Bin-discretized Y data as fake X
% % Get bin sizes equivalent to desired max and min noise level
% minbin = ceil(range(X,'all') / 5.9);
% maxbin = round(range(X, 'all') / 1);
% nbin = round(linspace(minbin, maxbin, n));
% % Preallocate
% MI_disc = cell(1,n);
% fakeX = cell(1,n);
% bin_edges = cell(1,n);
% bin = cell(1,n);
% MI_disc(:) = {zeros(length(noise), repeats)};
% fakeX(:) = {nan(size(X))};
% % Loop over number of bins to discretize by
% for i = 1:n
%     [~, bin_edges{i}, bin{i}] = histcounts(X, nbin(i));
%     fakeX{i} = bin_edges{i}(bin{i} + 1);
%     fakeX{i}(fakeX{i}==bin_edges{i}(1)) = nan;
%     % Indexing above can sometimes cause transpose, so undo if it happens
%     if size(fakeX{i}, 1) < size(fakeX{i}, 2) 
%         fakeX{i} = fakeX{i}';
%     end
%     MI_disc{i} = KSG_precision(fakeX{i}, Y, knn, repeats, noise);
% end
MI_real = KSG_precision(X, Y, knn, repeats, noise);

%% Get precision values
precision = zeros(1, n);
precision_ind = zeros(1, n);
for i = 1:n
    mis = MI_KSG_subsampling_multispike(fakeX{i}, Y, knn, (1:5));
    mi_sd = findMI_KSG_stddev(mis, size(X,1), false);
    precision_ind(i) = find(mean(MI_disc{i}, 2) < ((MI_disc{i}(1,1) - mi_sd)), 1);
    precision(i) = noise(precision_ind(i));
end

% Alternative ways to find precision
% Derivative threshold
deriv_thresh = -1.5e-3;
deriv_precision = zeros(1, n);
deriv_prec_ind = zeros(1, n);
for i = 1:n
    meanMI = mean(MI_disc{i}, 2);
    % Get smooth derivative of MI 
    [b,g] = sgolay(2, 15);
    % Note that magnitude of derivative is arbitrary
    grad = conv(meanMI, -1 * g(:,2), 'same'); 
    % Get where derivative passes below threshold
    deriv_prec_ind(i) = find(grad < deriv_thresh, 1);
    deriv_precision(i) = noise(deriv_prec_ind(i));
%     % Fit both lines, find intersection
%     b_start = [ones(start_ind, 1), noise(1:start_ind)'] \ meanMI(1:start_ind);
%     b_end = [ones(length(noise)-end_ind+1, 1), noise(end_ind:end)'] \ meanMI(end_ind:end);
end

%% Figures

% MI against noise at different binning levels
figure
hold on
cols = copper(n);
for i = 1:n
    % Plot line of MI against noise
%     mseb(log10(noise), ...
%         mean(MI_disc{i},2), ...
%         std(MI_disc{i}, 0, 2)',...
%         struct('col', {{cols(i,:)}}));
    plot(log10(noise), mean(MI_disc{i}, 2), 'color', cols(i,:))
    % Plot precision points
%     plot(log10(noise(precision_ind(i))), mean(MI_disc{i}(precision_ind(i),:)), '*')
%     plot(log10(noise(deriv_prec_ind(i))), mean(MI_disc{i}(deriv_prec_ind(i),:)), 'o')
end
mseb(log10(noise), mean(MI_real,2), std(MI_real, 0, 2)');
xlabel('log10(noise)')
ylabel('MI (bits)')
title(['Moth ', moth, ' ', muscle])

% 
figure
hold on
meanMI = mean(MI_real, 2);
for i = 1:n
%     [~, ind] = min(abs(mean(MI_disc{i}(1,:)) - mean(MI_real,2)));
%     plot(noise(ind), diff(bin_edges{i}(1:2)), '*')
    [~, ind] = min(abs(noise-diff(bin_edges{i}(1:2))));
    plot(meanMI(ind), MI_disc{i}(1,1), '*');
end
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
xlabel('MI at noise level equal to bin size')
ylabel('Initial MI at associated bin size')
title(['Moth ', moth, ' ', muscle])
%%
% Initial MI against number of bins
figure
hold on
entropy = @(n) - n .* (1 ./ n .* log2(1 ./ n));
plot(nbin, entropy(nbin) + (mean(MI_disc{1}(2,:)) - entropy(nbin(1))))
plot(nbin, cellfun(@(x) mean(x(1,1)), MI_disc), '.')
xlabel('# of bins')
ylabel('Initial MI (bits)')
title(['Moth ', moth, ' ', muscle])

% %%
% % Initial information vs prenoise (aka a priori precision)
% figure
% hold on
% mseb(log10(noise), mean(MI_prenoise{1},2), std(MI_prenoise{1}, 0, 2)');
% plot(log10(prenoise), cellfun(@(x) mean(x(1,:)), MI_prenoise), '*')
% 
% % Observed precision vs prenoise
% figure 
% hold on
% plot(prenoise, precision - precision(1), '*')
% plot(prenoise, deriv_precision - deriv_precision(1), 'o')
% %     plot(get(gca, 'xlim'), get(gca, 'xlim'), 'k-') 
% xlabel('Pre-added noise amplitude')
% ylabel('$\Delta$Precision observed', 'interpreter', 'latex')

%%

figure 
hold on
bob = 0.5 * rand(1000);
bob1 = 1 * rand(1000) + 0.5 * rand(1000);
bob2 = 1.5 * rand(1000);
histogram(bob, 'normalization', 'pdf')
histogram(bob1, 'normalization', 'pdf')
histogram(bob2, 'normalization', 'pdf')

%% Distribution of kth nearest neighbor distances as noise is added
% Levels of noise to check
ncheck = 5;
noiseinds = linspace(1, length(noise), ncheck);
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
        % Loop over noise levels
        for ii = 1:ncheck
            kdist{ii,jj} = nan(size(useX,1),1);
            noiseX = useX + noise(noiseinds(ii)) * rand(size(useX));
            noiseY = useY;
            % Z-score X and Y by column (same as in kraskov C code)
            xme = mean(noiseX, 1); 
            xsd = std(noiseX, 1);
            yme = mean(noiseY, 1); 
            ysd = std(noiseY, 1);
            noiseX = (noiseX - xme) .* xsd;
            noiseY = (noiseY - yme) .* ysd;
            % Loop over each point in X, get knn distances
            % (Can improve by passing all points to knnsearch in single call)
            for i = 1:size(useX,1)
                % Get indices of kth nearest neighbors, correct indices (X and Y)
                [idx, distx] = knnsearch(noiseX([1:(i-1), (i+1):end], :), noiseX(i,:), 'K', knn);
                idx(idx>i) = idx(idx>i) + 1;
                [idy, disty] = knnsearch(useY([1:(i-1), (i+1):end], :), useY(i,:), 'K', knn);
                idy(idy>i) = idy(idy>i) + 1;
                % Select max distance, put back in regular units
                if distx(knn) > disty(knn)
                    kdist{ii,jj}(i) = sum((((noiseX(idx(4),:)/xsd+xme) - (noiseX(i,:)/xsd+xme))).^2).^0.5;
                else
                    kdist{ii,jj}(i) = sum((((noiseY(idy(4),:)/ysd+yme) - (noiseY(i,:)/ysd+yme))).^2).^0.5;
                end
            end
            % Make plot
            ax(ii) = subaxis(ncheck, 1 , ii, 'SpacingVert', 0);
            histogram(log10(kdist{ii,jj}), 40, 'Normalization', 'pdf', 'EdgeColor', 'none')
            xline(median(log10(kdist{ii,jj}), 'omitnan'), 'k')
            text(0.6, 0.5, ['r = ', num2str(noise(noiseinds(ii)))], 'units', 'normalized')
            xlims(ii,:) = get(gca, 'Xlim');
        end
        % Link axes so X axis is shared, clean up X ticks
        xtickspace = logspace(min(xlims(:,1)), max(xlims(:,2)), 8);
        for ii = 1:ncheck
            ax(ii).set('xlim', [min(xlims(:,1)), max(xlims(:,2))]);
            ax(ii).set('xtick', log10(xtickspace));
            ax(ii).set('xticklabels', []);
        end
        ax(end).set('xticklabels', num2cell(round(10.^get(ax(end),'XTick'), 4)));
        xlabel(ax(end), ...
            ['k=',num2str(knn),'$^{th}$ N.N. distance in space $\|z-z`\| = max\{\|x-x`\|, \|y-y`\|\}$'], ...
            'interpreter', 'latex')
        ylabel(ax(round(length(ax)/2)), 'Probability')
        title(ax(1), ['Moth ', moth,', ', muscle, ', ', num2str(jj), ' spikes in wingbeat'])
    end
end

% Median distance against noise
medians = cellfun(@(x) median(x, 'omitnan'), kdist);
figure()
plot(noise(noiseinds), medians, '.-', 'markersize', 10)
xlabel('$r_c$ (ms)', 'interpreter', 'latex')
ylabel('Median knn distance (ms)')
title(['Moth ', moth,', ', muscle])
enough_data = arrayfun(@(x) sum(Nspike==x)>knn, unq);
legend(arrayfun(@(x) [num2str(x), ' spike(s)'], unq(unq~=0 & enough_data), 'UniformOutput', false), ...
    'location', 'northwest')

% Take a rough guess at what the overall scale at zero noise is
probs = arrayfun(@(x) sum(Nspike==x)/length(Nspike), unq);
probs(unq==0) = [];
scale = sum(probs' .* medians(1,:), 'omitnan') / sum(probs,'omitnan')
