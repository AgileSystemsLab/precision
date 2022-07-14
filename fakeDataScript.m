rng('shuffle')
% Moth and muscle to focus on
moth = '3';
muscle = 'LDVM';
knn = 4;
noise = (0:.05:6);
repeats = 150;
n = 6; % How many different levels of rescaling to use

%---- Load data
load(fullfile('Data',['Moth',num2str(moth),'_MIdata.mat']))
% Rename so it's easier to write, get some useful quantities out
X = time_data.([muscle,'strokes']);
Y = Tz_WSd;
Nspike = sum(~isnan(X), 2);
unq = unique(Nspike);

%---- Run precision estimation on real data
% MI_real = KSG_precision(X, Y, knn, repeats, noise, true);
% title(['Real data, moth ', moth, ', ', muscle])

% %---- Fake data: Totally uncorrelated
% % Same resolution, range, and # of spikes as real data, but uniformly distributed
% fake_uncorr = nan(size(X));
% % Shift so min resolution is 1, allowing rand generation using ints
% minmax = round([min(X, [], 'all'), max(X, [], 'all')] * 10); 
% for i = 1:size(X,1)
%     fake_uncorr(i,1:Nspike(i)) = randi(minmax, 1, Nspike(i));
% end
% fake_uncorr = fake_uncorr / 10;

%---- Direct connected fake data ahead-of-time corrupted to specific noise levels
prenoise = linspace(0, 6, 5);
n = 5;
% fakeY = normrnd(0, 1, 3000, 2);
% fakeX = fakeY * range(X,'all') / range(fakeY, 'all');
[MI_prenoise, precision, precision_ind] = precorruption_analysis(X, Y, 4, n, prenoise, noise, 150);

%%
% Precision curves
figure
hold on
cols = copper(n);
for i = 1:n
    rescale = 1 / mean(MI_prenoise{i}(1,:));
%     rescale = 1;
    mseb(log10(noise), ...
        rescale * mean(MI_prenoise{i},2), rescale * std(MI_prenoise{i}, 0, 2)',...
        struct('col', {{cols(i,:)}}));
    plot(log10(noise(precision_ind(i))), mean(MI_prenoise{i}(precision_ind(i),:)), '*')
end

% Initial information vs prenoise (aka a priori precision)
figure
hold on
mseb(log10(noise), mean(MI_prenoise{1},2), std(MI_prenoise{1}, 0, 2)');
plot(log10(prenoise), cellfun(@(x) mean(x(1,:)), MI_prenoise), '*')

% Observed precision vs prenoise
figure 
hold on
plot(precision(1) + prenoise, precision - precision(1), '*')
%     plot(get(gca, 'xlim'), get(gca, 'xlim'), 'k-')
xlabel('Pre-added noise amplitude')
ylabel('$\Delta$Precision observed', 'interpreter', 'latex')


% %---- Deterministically linked fake data at different sample sizes
% figure()
% hold on
% cols = copper(n);
% for i = 1:n
%     signal_amp = max(fake_X{i}, [], 'all') - min(fake_X{i}, [], 'all');
%     rescale = 1 / mean(MI_direct{i}(1,:));
% %     rescale = 1;
%     plot(log10(noise), rescale * mean(MI_direct{i}, 2), 'color', cols(i,:));
% end
% title('motor PCs directly to spike times with varying sample size')
% xlabel('log10( noise added / signal range )')
% ylabel('normalized information (1 at zero noise added)')
% colormap(copper)
% cbh = colorbar;
% set(cbh,'YTick',linspace(0,1,n))
% set(cbh,'YTickLabel', num2str(data_fraction.'))



% %% Distribution of kth nearest neighbor distances as noise is added
% % Levels of noise to check
% ncheck = 5;
% kdist = cell(ncheck, length(unq(unq~=0)));
% % Loop over number of spikes in a wingbeat
% % for jj = unq(unq~=0)'
% for jj = 2
% %     useX = X(Nspike==jj, 1:jj);
% %     useY = Y(Nspike==jj, :);
%     useX = Y;
%     useY = Y;
%     % Continue only if enough wingbeats with this many spikes
%     if size(useX, 1) >= knn
%         % Preallocate, setup figure
%         figure
%         hold on
%         xlims = zeros(ncheck, 2);
%         ax = gobjects(ncheck, 1);
%         noiseinds = linspace(1, length(noise), ncheck);
%         % Loop over noise levels
%         for ii = 1:ncheck
%             noiseamp = noise(noiseinds(ii));
%             kdist{ii,jj} = nan(size(useX,1),1);
%             noiseX = useX + noiseamp * rand(size(useX));
%             noiseY = useY;
%             % Z-score X and Y by column (same as in kraskov C code)
%             xme = mean(noiseX, 1); 
%             xsd = std(noiseX, 1);
%             yme = mean(noiseY, 1); 
%             ysd = std(noiseY, 1);
%             noiseX = (noiseX - xme) .* xsd;
%             noiseY = (noiseY - yme) .* ysd;
%             % Loop over each point in X, get knn distances
%             % (Can improve by passing all points to knnsearch in single call)
%             for i = 1:size(useX,1)
%                 % Get indices of kth nearest neighbors, correct indices (X and Y)
%                 [idx, distx] = knnsearch(noiseX([1:(i-1), (i+1):end], :), noiseX(i,:), 'K', knn);
%                 idx(idx>i) = idx(idx>i) + 1;
%                 [idy, disty] = knnsearch(useY([1:(i-1), (i+1):end], :), useY(i,:), 'K', knn);
%                 idy(idy>i) = idy(idy>i) + 1;
%                 % Select max distance
% %                 kdist{ii,jj}(i) = max(distx(knn), disty(knn));
%                 % Put back in regular units
%                 if distx(knn) > disty(knn)
%                     kdist{ii,jj}(i) = sum((((noiseX(idx(4),:)/xsd+xme) - (noiseX(i,:)/xsd+xme))).^2).^0.5;
%                 else
%                     kdist{ii,jj}(i) = sum((((noiseY(idy(4),:)/ysd+yme) - (noiseY(i,:)/ysd+yme))).^2).^0.5;
%                 end
%             end
%             % Make plot
%             ax(ii) = subaxis(ncheck, 1 , ii, 'SpacingVert', 0);
%             histogram(log10(kdist{ii,jj}), 40, 'Normalization', 'pdf', 'EdgeColor', 'none')
%             xline(median(log10(kdist{ii,jj}), 'omitnan'), 'k')
%             text(0.6, 0.5, ['r = ', num2str(noiseamp)], 'units', 'normalized')
%             xlims(ii,:) = get(gca, 'Xlim');
%         end
%         % Link axes so X axis is shared, clean up X ticks
% %         xtickspace = linspace(10^min(xlims(:,1)), 10^max(xlims(:,2)), 8);
%         xtickspace = logspace(min(xlims(:,1)), max(xlims(:,2)), 8);
%         for ii = 1:ncheck
%             ax(ii).set('xlim', [min(xlims(:,1)), max(xlims(:,2))]);
%             ax(ii).set('xtick', log10(xtickspace));
%             ax(ii).set('xticklabels', []);
%         end
%         ax(end).set('xticklabels', num2cell(round(10.^get(ax(end),'XTick'), 2)));
%         xlabel(ax(end), ...
%             ['k=',num2str(knn),'$^{th}$ N.N. distance in space $\|z-z`\| = max\{\|x-x`\|, \|y-y`\|\}$'], ...
%             "interpreter", "latex")
%         ylabel(ax(round(length(ax)/2)), 'Probability')
%         title(ax(1), ['Moth ', moth,', ', muscle, ', ', num2str(jj), ' spikes in wingbeat'])
%     end
% end
% 
% %% Median distance against noise
% medians = cellfun(@(x) median(x, 'omitnan'), kdist);
% figure()
% boxplot(medians)
% xlabel('# of spikes in wingbeat')
% ylabel('Median knn distance across noise levels (ms)')
% title(['Moth ', moth,', ', muscle])
% 
% % Take a rough guess at what the overall scale at zero noise is
% probs = arrayfun(@(x) sum(Nspike==x)/length(Nspike), unq);
% probs(unq==0) = [];
% scale = sum(probs' .* medians(1,:), 'omitnan') / sum(probs,'omitnan')


function [MI_prenoise, precision, precision_ind] = precorruption_analysis(X, Y, knn, n, prenoise, noise, repeats)
    arguments
        X (:,:) double
        Y (:,:) double
        knn (1,1) double = 4
        n (1,1) double = 5
        prenoise (1,:) double = linspace(0, 6, n)
        noise (1,:) double = (0:0.05:6)
        repeats (1,1) double = 100
    end
    % Setup parameters and preallocate
    MI_prenoise = cell(1,n);
    MI_prenoise(:) = {zeros(length(noise), 2)};
    precision = zeros(1, n);
    precision_ind = zeros(1, n);
    % Loop over noise levels, create fake data, run precision estimation
    for i = 1:n
        fake_X = X + prenoise(i) * rand(size(X));
        MI_prenoise{i} = KSG_precision(fake_X, Y, knn, repeats, noise);
        % Get precision value
        mis = MI_KSG_subsampling_multispike(fake_X, Y, knn, (1:50));
        mi_sd = findMI_KSG_stddev(mis, size(X,1), false);
        precision_ind(i) = find(mean(MI_prenoise{i}, 2) < ((mean(MI_prenoise{i}(1,:)) - 2*mi_sd)), 1);
        precision(i) = noise(precision_ind(i));
    end
    
%     % Precision curves
%     figure
%     hold on
%     cols = copper(n);
%     for i = 1:n
%     %     rescale = 1 / mean(MI_prenoise{i}(1,:));
%         rescale = 1;
%         mseb(log10(noise), ...
%             rescale * mean(MI_prenoise{i},2), rescale * std(MI_prenoise{i}, 0, 2)',...
%             struct('col', {{cols(i,:)}}));
%         plot(log10(noise(precision_ind(i))), mean(MI_prenoise{i}(precision_ind(i),:)), '*')
%     end
%     
%     % Initial information vs prenoise (aka a priori precision)
%     figure
%     hold on
%     mseb(log10(noise), mean(MI_prenoise{1},2), std(MI_prenoise{1}, 0, 2)');
%     plot(log10(prenoise), cellfun(@(x) mean(x(1,:)), MI_prenoise), '*')
%     
%     % Observed precision vs prenoise
%     figure 
%     hold on
%     plot(precision(1) + prenoise, precision - precision(1), '*')
% %     plot(get(gca, 'xlim'), get(gca, 'xlim'), 'k-')
%     xlabel('Pre-added noise amplitude')
%     ylabel('$\Delta$Precision observed', 'interpreter', 'latex')
end