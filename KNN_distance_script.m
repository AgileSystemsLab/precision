% Distribution of kth nearest neighbor distances as noise is added

run_all = false;


%% Detailed plots of single example
% Moth and muscle to use as example
moth = '3';
muscle = 'RSA';
% Note: Moth 3 RSA precision with STD method is 0.5368 ms
knn = 4;
noise = [0, logspace(log10(0.05), log10(6), 120)];
repeats = 150;
n = 10; 

%---- Load data
load(fullfile('Data',['Moth',num2str(moth),'_MIdata.mat']))
fields = fieldnames(time_data);
% Rename so it's easier to write, get some useful quantities out
X = time_data.([muscle,'strokes']);
Y = Tz_WSd;
Nspike = sum(~isnan(X), 2);
unq = unique(Nspike);

% Figure setup
figure
hold on
cols = {'#e41a1c', '#377eb8', '#4daf4a'};
% Levels of noise to check
noiselevels = [0, 0.75, 1.5, 3, 6];
ncheck = length(noiselevels);
kdist = cell(ncheck, length(unq(unq~=0)));
% Loop over number of spikes in a wingbeat
for jj = unq(unq~=0)'
    useX = X(Nspike==jj, 1:jj);
    useY = Y(Nspike==jj, :);
    % Continue only if enough wingbeats with this many spikes
    if size(useX, 1) >= knn
        % Preallocate, setup figure
        xlims = zeros(ncheck, 2);
        ax = gobjects(ncheck, 1);
        % Loop over noise levels
        for ii = 1:ncheck
            kdist{ii,jj} = nan(size(useX,1),1);
            noiseX = useX + noiselevels(ii) * rand(size(useX));
            noiseY = useY;
            % Add noise just as in kraskov C code
            noiseX = noiseX + 1e-8*rand(size(noiseX));
            noiseY = noiseY + 1e-8*rand(size(noiseY));
            % Z-score X and Y by column (same as in kraskov C code)
            xme = mean(noiseX, 1); 
            xsd = std(noiseX, 1);
            yme = mean(noiseY, 1); 
            ysd = std(noiseY, 1);
            noiseX = (noiseX - xme) ./ xsd;
            noiseY = (noiseY - yme) ./ ysd;
            % Loop over each point in X, get knn distances
            % (Can improve by passing all points to knnsearch in single call)
            for i = 1:size(useX,1)
                % Get indices of kth nearest neighbors, correct indices (X and Y)
                [idx, ~] = knnsearch(noiseX([1:(i-1), (i+1):end], :), noiseX(i,:), 'K', knn, 'NSMethod', 'kdtree');
                idx(idx>i) = idx(idx>i) + 1;
                [idy, ~] = knnsearch(useY([1:(i-1), (i+1):end], :), useY(i,:), 'K', knn, 'NSMethod', 'kdtree');
                idy(idy>i) = idy(idy>i) + 1;
                % Calculate max norm distances
                distx = max(sqrt((noiseX(i,:) - noiseX(idx(knn),:)).^2));
                disty = max(sqrt((noiseY(i,:) - noiseY(idy(knn),:)).^2));
                % Select max distance, put back in regular units
                if distx > disty
                    kdist{ii,jj}(i) = distx;
                else
                    kdist{ii,jj}(i) = disty;
                end
            end
            % Make plot
            ax(ii) = subaxis(ncheck, 1 , ii, 'SpacingVert', 0);
            hold on
            binedges = logspace(log10(0.01), log10(3.5), 40);
            histogram(kdist{ii,jj}, binedges, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', cols{jj}, 'FaceAlpha', 0.5)
            xline(median(kdist{ii,jj}, 'omitnan'), 'color', cols{jj}, 'LineWidth', 1)
            text(0.7, 0.5, ['r = ', num2str(noiselevels(ii))], 'units', 'normalized')
            xlims(ii,:) = get(gca, 'Xlim');
        end
        % Link axes so X axis is shared, clean up X ticks
        xtickspace = logspace(-2, log10(max(xlims(:,2))), 5);
        for ii = 1:ncheck
            ax(ii).set('Xscale', 'log')
            ax(ii).set('xlim', [min(xlims(:,1)), max(xlims(:,2))]);
            ax(ii).set('xtick', xtickspace);
            if ii ~= ncheck
                ax(ii).set('xticklabels', []);
            end
        end
        xlabel(ax(end), ...
            ['k=',num2str(knn),'$^{th}$ N.N. distance in space $\|z-z`\| = max\{\|x-x`\|, \|y-y`\|\}$ (a.u.)'], ...
            'interpreter', 'latex')
        ylabel(ax(round(length(ax)/2)), 'Probability')
        title(ax(1), ['Moth ', moth,' ', muscle])
        legend(ax(3), {'1 spike', '', '2 spikes', '', '3 spikes', ''})
    end
end
exportgraphics(gcf, fullfile('figures','KNN_distances_example_distributions.pdf'),'ContentType','vector')

% Median distance against noise
medians = cellfun(@(x) median(x, 'omitnan'), kdist);
figure('OuterPosition', [962, 377, 385, 499])
hold on
for i = 1:size(medians, 2)
    plot(noiselevels, medians(:,i), '.-', 'markersize', 10, 'color', cols{i})
end
xlabel('$r_c$ (ms)', 'interpreter', 'latex')
ylabel('Median knn distance (a.u.)')
title(['Moth ', moth,', ', muscle])
enough_data = arrayfun(@(x) sum(Nspike==x)>knn, unq);
legend(arrayfun(@(x) [num2str(x), ' spike(s)'], unq(unq~=0 & enough_data), 'UniformOutput', false), ...
    'location', 'northwest')
exportgraphics(gcf, fullfile('figures','KNN_distances_example_medians.pdf'),'ContentType','vector')



%% Run ALL
if run_all
    nmoths = 7;
    nmuscles = 10;
    
    % Levels of noise to check
    noiselevels = [0, 0.75, 1.5, 3, 6];
    ncheck = length(noiselevels);
    
    % Preallocate
    kdist = cell(nmoths, nmuscles);
    xme = cell(nmoths, nmuscles);
    yme = cell(nmoths, nmuscles);
    xsd = cell(nmoths, nmuscles);
    ysd = cell(nmoths, nmuscles);
    
    for mothi = 1:nmoths
        %---- Load data
        load(fullfile('Data',['Moth',num2str(mothi),'_MIdata.mat']))
        fields = fieldnames(time_data);
        disp(['Moth ', num2str(mothi)])
        for musci = 1:nmuscles
            disp(fields{musci}(1:end-7))
            tic
            % Rename so it's easier to write, get some useful quantities out
            X = time_data.(fields{musci});
            Y = Tz_WSd;
            Nspike = sum(~isnan(X), 2);
            unq = unique(Nspike);
            % Preallocate
            kdist{mothi,musci} = cell(ncheck, length(unq(unq~=0)));
            xme{mothi,musci} = cell(ncheck, length(unq(unq~=0)));
            yme{mothi,musci} = cell(ncheck, length(unq(unq~=0)));
            xsd{mothi,musci} = cell(ncheck, length(unq(unq~=0)));
            ysd{mothi,musci} = cell(ncheck, length(unq(unq~=0)));
            % Loop over number of spikes in a wingbeat
            for jj = unq(unq~=0)'
                useX = X(Nspike==jj, 1:jj);
                useY = Y(Nspike==jj, :);
                % Continue only if enough wingbeats (knn+1 or more) with this many spikes
                if size(useX, 1) > knn
                    % Loop over noise levels
                    for ii = 1:ncheck
                        kdist{mothi,musci}{ii,jj} = nan(size(useX,1),1);
                        noiseX = useX + noiselevels(ii) * rand(size(useX));
                        noiseY = useY;
                        % Z-score X and Y by column (same as in kraskov C code)
                        xme{mothi,musci}{ii,jj} = mean(noiseX, 1); 
                        xsd{mothi,musci}{ii,jj} = std(noiseX, 1);
                        yme{mothi,musci}{ii,jj} = mean(noiseY, 1); 
                        ysd{mothi,musci}{ii,jj} = std(noiseY, 1);
                        noiseX = (noiseX - xme{mothi,musci}{ii,jj}) ./ xsd{mothi,musci}{ii,jj};
                        noiseY = (noiseY - yme{mothi,musci}{ii,jj}) ./ ysd{mothi,musci}{ii,jj};
                        % Loop over each point in X, get knn distances
                        % (Can improve by passing all points to knnsearch in single call)
                        for i = 1:size(useX,1)
                            % Get indices of kth nearest neighbors, correct indices (X and Y)
                            [idx, ~] = knnsearch(noiseX([1:(i-1), (i+1):end], :), noiseX(i,:), 'K', knn, 'NSMethod', 'kdtree');
                            idx(idx>i) = idx(idx>i) + 1;
                            [idy, ~] = knnsearch(useY([1:(i-1), (i+1):end], :), useY(i,:), 'K', knn, 'NSMethod', 'kdtree');
                            idy(idy>i) = idy(idy>i) + 1;
                            % Calculate max norm distances
                            distx = max(sqrt((noiseX(i,:) - noiseX(idx(knn),:)).^2));
                            disty = max(sqrt((noiseY(i,:) - noiseY(idy(knn),:)).^2));
                            % Select max distance, leave in z-scored units
                            if distx > disty
                                kdist{mothi,musci}{ii,jj}(i) = distx;
                            else
                                kdist{mothi,musci}{ii,jj}(i) = disty;
                            end
                        end
                    end
                end
            end
            toc
        end
    end
    
    save('KSG_data_knn_distances.mat', 'kdist', 'xme', 'yme', 'xsd', 'ysd', 'noiselevels', 'ncheck');
end

%% Plot results of all distances
load('KSG_data_knn_distances.mat')

nmoths = 5;
nmuscles = 10;

%--- Change in median knn distance
allmedians = cell(7,1);
for i = 1:nmoths
    for j = 1:nmuscles
        un_zscored = cellfun(@(x,y) x * y, ...
            kdist{i,j}, ...
            cellfun(@(x) mean(x), xsd{i,j}, 'UniformOutput', false), ...
            'UniformOutput', false);
        medians = cellfun(@(x) median(x, 'omitnan'), un_zscored);
        % Loop over number of spikes
        for k = 1:size(medians, 2)
            allmedians{k} = [allmedians{k}; medians(:,k)' - medians(1,k)];
        end
    end
end
figure('OuterPosition', [590 258 900 500])
hold on
grid on
cols = copper(7);
for i = 1:7
    if size(allmedians{i}, 1) > 1
        errorbar(noiselevels-i*0.08+0.24, mean(allmedians{i}, 1, 'omitnan'), std(allmedians{i}, 1, 'omitnan'), ...
            'color', cols(i,:), 'Marker', 'o', 'MarkerEdgeColor', cols(i,:), 'MarkerFaceColor', cols(i,:), ...
            'LineWidth', 1, ...
            'DisplayName', [num2str(i), ' spikes'])
    else
        plot(noiselevels, allmedians{i}, '.-', 'color', cols(i,:), ...
            'Marker', 'o', 'MarkerEdgeColor', cols(i,:), 'MarkerFaceColor', 'w', ...
            'LineWidth', 1, 'LineStyle', '--', ...
            'DisplayName', [num2str(i), ' spikes'])
    end
end
set(gca, 'XTick', noiselevels)
legend('Location', 'northwest')
ylabel('Change in median k-NN distance from zero noise (ms)')
xlabel('Noise amplitude r_c (ms)')
exportgraphics(gcf, fullfile('figures','KNN_distances_all_medians.pdf'),'ContentType','vector')


figure('OuterPosition', [724 77 590 783])
cols = copper(7);
ax = gobjects(5,2);
for i = 1:5
    for j = 1:2
        ax(i,j) = subaxis(5, 2, j, i, 'SpacingVert', 0.01, 'SpacingHoriz', 0.05);
        hold on
    end
end

for j = 1:nmuscles
    allmedians = cell(7,1);
    for i = 1:nmoths
        un_zscored = cellfun(@(x,y) x * y, ...
            kdist{i,j}, ...
            cellfun(@(x) mean(x), xsd{i,j}, 'UniformOutput', false), ...
            'UniformOutput', false);
        medians = cellfun(@(x) median(x, 'omitnan'), un_zscored);
        % Loop over number of spikes
        for k = 1:size(medians, 2)
            allmedians{k} = [allmedians{k}; medians(:,k)' - medians(1,k)];
        end
    end
    row = mod(j-1,5)+1;
    col = int8(j>5)+1;
    set(gcf, 'CurrentAxes', ax(row,col))
    for i = 1:7
        if size(allmedians{i}, 1) > 1
            errorbar(noiselevels-i*0.08+0.24, mean(allmedians{i}, 1, 'omitnan'), std(allmedians{i}, 1, 'omitnan'), ...
                'color', cols(i,:), 'Marker', 'o', 'MarkerEdgeColor', cols(i,:), 'MarkerFaceColor', cols(i,:), ...
                'LineWidth', 1, ...
                'DisplayName', [num2str(i), ' spikes'])
        elseif ~isempty(allmedians{i})
            plot(noiselevels, allmedians{i}, '.-', 'color', cols(i,:), ...
                'Marker', 'o', 'MarkerEdgeColor', cols(i,:), 'MarkerFaceColor', 'w', ...
                'LineWidth', 1, 'LineStyle', '--', ...
                'DisplayName', [num2str(i), ' spikes'])
        end
    end
    text(0.1, 0.9, fields{j}(1:end-7), 'units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold')
    ylim([-0.1, 2.1])
    set(gca, 'XTick', noiselevels)
end
% Global axes (huge pain, WHY MATLAB)
p1 = get(ax(1,1), 'position');
p2 = get(ax(end,1), 'position');
p3 = get(ax(end,end), 'position');
height = p1(2) + p1(4) - p2(2);
width = p3(1) + p3(3) - p2(1);
axh = axes('position', [p2(1), p2(2), width, height], 'visible', 'off');
axh.XLabel.Visible = 'on';
axh.YLabel.Visible = 'on';
axes(axh)
xlabel('Noise amplitude r_c (ms)')
ylabel('Change in median k-NN distance from zero noise (ms)')
exportgraphics(gcf, fullfile('figures','KNN_distances_all_medians_by_muscle.pdf'),'ContentType','vector')