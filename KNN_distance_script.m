% Distribution of kth nearest neighbor distances as noise is added
% Moth and muscle to focus on
moth = '3';
muscle = 'RSA';
knn = 4;
noise = [0, logspace(log10(0.05), log10(6), 120)];
repeats = 150;
n = 10; 

%---- Load data
load(fullfile('Data',['Moth',num2str(moth),'_MIdata.mat']))
% Rename so it's easier to write, get some useful quantities out
X = time_data.([muscle,'strokes']);
Y = Tz_WSd;
Nspike = sum(~isnan(X), 2);
unq = unique(Nspike);

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
                    kdist{ii,jj}(i) = sum((((noiseX(idx(knn),:)/xsd+xme) - (noiseX(i,:)/xsd+xme))).^2).^0.5;
                else
                    kdist{ii,jj}(i) = sum((((noiseY(idy(knn),:)/ysd+yme) - (noiseY(i,:)/ysd+yme))).^2).^0.5;
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