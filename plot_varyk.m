% Plot precision while varying k

load('KSG_data_varyk.mat')

nmoths = 7;
nmuscles = 10;
nk = length(knn_vec);
nsubsets = 6;

precision = nan(nmoths, nmuscles, nk);
for i = 1:nmoths
    % Load data (for finding precision values)
    load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
    fields = fieldnames(time_data);
    for j = 1:nmuscles
        for k = 1:nk
            knn = knn_vec(k);
            % Get precision and plot
            meanMI = mean(MI_varyk{i,j,k}, 2);
            mis = MI_KSG_subsampling_multispike(time_data.(fields{j}), Tz_WSd, knn, (1:nsubsets));
            mi_sd = findMI_KSG_stddev(mis, size(Tz_WSd,1), false);
            ind = find(meanMI < (MI_varyk{i,j,k}(1,1) - mi_sd), 1);
            if ~isempty(ind)
                precision(i,j,k) = noise(ind);
            end
        end
    end
end

%%
figure('OuterPosition', [669 141 403 677])
ax = gobjects(5,1);
col = '#4472C4';

for i = 1:5
    ax(i) = subaxis(5, 1 , i, 'SpacingVert', 0.02);
    hold on
    pr = [reshape(precision(:,i,:), nmoths, nk); reshape(precision(:,5+i,:), nmoths, nk)];
    errorbar(knn_vec, mean(pr, 1), std(pr, 1), 'o', ...
        'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
        'LineWidth', 1, 'CapSize', 13)
    for j = 1:nk
        plot(knn_vec(j) + 0.1 * normrnd(0, 0.5, 1, 2*nmoths), pr(:,j), 'k.')
    end
    text(0.1, 0.9, fields{i}(2:end-7), 'units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold')
    xl = get(gca, 'xlim');
    xlim([xl(1)-0.5, xl(2)+0.5])
    ylim([0 4])
    set(ax(i), 'XTick', knn_vec)
    set(ax(i), 'Xticklabel', [])
    if i == 3
        ylabel('Spike Timing Precision r_c (ms)')
    elseif i == 5
        xlabel('K')
        set(ax(i), 'Xticklabel', knn_vec)
    end
end

% Save figure
exportgraphics(gcf,fullfile('figures','KSG_precision_vary_with_K.pdf'),'ContentType','vector')