% Plot all KSG and NSB

% Load data
load('KSG_data_count.mat')
load('KSG_data.mat')
load('NSB_data.mat')
load('NSB_bias_data.mat')

nmoths = 7;
nmuscles = 10;
% NSB parameters
ntorque = 3;
nspikebins = 70;

% STD method
nsubsets = 8; 

muscle_names = {'L3AX'; 'LBA'; 'LSA'; 'LDVM'; 'LDLM'; 'R3AX'; 'RBA'; 'RSA'; 'RDVM'; 'RDLM'};

%%
load('NSB_data.mat')
load('NSB_bias_data.mat')

STD = sqrt(dS_nsbwordvec.^2 + conditionaldS_nsbvec(:, :, ntorque).^2);
bias = mean(conditionalentropyvec_bias, 4);

% Set up figure
figure
hold on
cols = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
% Loop over moths
for i = 1:nmoths
    load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
    fields = fieldnames(time_data);
    j = 1;
    % Get bin sizes
    bins = range(time_data.(fields{j}), 'all') ./ (1:nspikebins);
    % Plot 
    ind = (i-1)*10 + j;
    % Plot main precision line
    y_bias = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
    y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque);
    mseb(bins, y, STD(ind, :), struct('col', {{cols{i}}}), 1);
    plot(bins, y_bias, '--', 'color', cols{i}, 'LineWidth', 3)
end
set(gca, 'Xscale', 'log')

%% Plot all NSB results
STD = sqrt(dS_nsbwordvec.^2 + conditionaldS_nsbvec(:, :, ntorque).^2);
bias = mean(conditionalentropyvec_bias, 4);

% Set up figure
figure('Outerposition', [597, 61, 833, 812])
ax = gobjects(5, 2);
for i = 1:5
    for j = 1:2
        ax(i,j) = subaxis(5, 2, j, i, 'SpacingVert', 0.01, 'SpacingHoriz', 0.03);
        hold on
    end
end
cols = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
precision_d = zeros(nmoths, nmuscles);
% Loop over moths
for i = 1:nmoths
    load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
    fields = fieldnames(time_data);
    % Loop over muscles
    for j = 1:nmuscles
        % Get bin sizes
        bins = range(time_data.(fields{j}), 'all') ./ (1:nspikebins);
        % Plot 
        row = mod(j-1,5)+1;
        col = int8(j>5)+1;
        ind = (i-1)*10 + j;
        % Plot main precision line
        y_bias = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
        y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque);
        set(gcf, 'CurrentAxes', ax(row,col))
%         mseb(bins, y, STD(ind, :), struct('col', {{cols{i}}}), 1);
        plot(bins, y_bias, 'color', cols{i}, 'LineWidth', 3)
        % Plot precision dots
        compvalue = mean(y_bias(end-10:end));
        compstd = std(y_bias(end-10:end));
        maxval = max(y_bias(y_bias>=compvalue));
        % Peak case, find max value location
        if maxval/compvalue >= 1.2
            [~,precision_ind] = max(y_bias);
        % Plateau case, find farthest right value near compvalue
        else
            precision_ind = find(y_bias >= (compvalue-2*compstd), 1);
        end
        plot(bins(precision_ind), y_bias(precision_ind), '.', 'MarkerSize', 30, 'Color', cols{i}, ...
            'DisplayName', 'no')
        plot(bins(precision_ind), y_bias(precision_ind), 'k.', 'MarkerSize', 15, 'DisplayName', 'no')
        precision_d(i,j) = bins(precision_ind);
    end
end
% Plot labels and aesthetics
for j = 1:nmuscles
    % Get row and column, set that ax to gca
    row = mod(j-1,5)+1;
    col = int8(j>5)+1;
    set(gcf, 'CurrentAxes', ax(row,col))
    % Subplot title
    text(0.1, 0.9, fields{j}(1:end-7), 'units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold')
    % Set all to log scale
    set(ax(row,col), 'Xscale', 'log')
    % Remove xticks on all but bottom
    if row~=5
        ax(row,col).set('xticklabels', []);
    end
    % Remove yticks on right column
    if col==2
        ax(row,col).set('yticklabels', []);
    end
    if (row==5) && (col==2)
        set_leg_off = findobj('DisplayName', 'no');
        for k = 1:numel(set_leg_off)
            set_leg_off(k).Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        lgd = legend(arrayfun(@(x) ['Moth ',num2str(x)], 1:nmoths, 'UniformOutput', false));
        lgd.BoxFace.ColorType = 'truecoloralpha';
        lgd.BoxFace.ColorData = uint8(255*[1 1 1 0.25]');
    end
    % Default all text on this axis to be larger
    set(ax(row,col), 'FontSize', 12)
end
% Set all to same x,y limits
linkaxes(ax)
% X and Y labels
han = axes(gcf,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('r_d (ms)', 'fontsize', 15)
han.XLabel.Position(2) = -0.05;
ylabel('MI (bits / wing stroke)', 'fontsize', 15)
han.YLabel.Position(1) = -0.08;
% Save figure
exportgraphics(gcf,fullfile('figures','NSB_all.pdf'),'ContentType','vector')

%% Plot all KSG results

% Set up figure
figure('Outerposition', [597, 61, 833, 812])
ax = gobjects(5, 2);
for i = 1:5
    for j = 1:2
        ax(i,j) = subaxis(5, 2, j, i, 'SpacingVert', 0.01, 'SpacingHoriz', 0.01);
        hold on
    end
end
cols = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
% Plot precision curves
for i = 1:nmoths
    for j = 1:nmuscles
        % Get row and column, set that ax to gca
        row = mod(j-1,5)+1;
        col = int8(j>5)+1;
        set(gcf, 'CurrentAxes', ax(row,col))
        % Plot!
        mseb(noise, mean(MI{i,j}, 2) + MI_count(i,j), std(MI{i,j}, 0, 2)', struct('col', {{cols{i}}}), 1);
    end
end
% Plot precision points
precision = nan(nmoths, nmuscles);
for i = 1:nmoths
    % Load data (for finding precision values)
    load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
    fields = fieldnames(time_data);
    for j = 1:nmuscles
        % Get row and column, set that ax to gca
        row = mod(j-1,5)+1;
        col = int8(j>5)+1;
        set(gcf, 'CurrentAxes', ax(row,col))
        % Get precision and plot
        meanMI = mean(MI{i,j}, 2);
        mis = MI_KSG_subsampling_multispike(time_data.(fields{j}), Tz_WSd, knn, (1:nsubsets));
        mi_sd = findMI_KSG_stddev(mis, size(Tz_WSd,1), false);
        ind = find(meanMI < (MI{i,j}(1,1) - mi_sd), 1);
        if isempty(ind)
            precision(i,j) = nan;
        else
            precision(i,j) = noise(ind);
            plot(noise(ind), meanMI(ind) + MI_count(i,j), ...
                '.', 'color', cols{i}, 'MarkerSize', 18)
            plot(noise(ind), meanMI(ind) + MI_count(i,j), ...
                'k.', 'MarkerSize', 10)
        end
    end
end
% Plot labels and aesthetics
for j = 1:nmuscles
    % Get row and column, set that ax to gca
    row = mod(j-1,5)+1;
    col = int8(j>5)+1;
    set(gcf, 'CurrentAxes', ax(row,col))
    % Subplot title
    text(0.1, 0.9, fields{j}(1:end-7), 'units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold')
    % Set all to log scale
    set(ax(row,col), 'Xscale', 'log')
    % Remove xticks on all but bottom
    if row~=5
        ax(row,col).set('xticklabels', []);
    end
    % Remove yticks on right column
    if col==2
        ax(row,col).set('yticklabels', []);
    end
    if (row==5) && (col==2)
        lgd = legend(arrayfun(@(x) ['Moth ',num2str(x)], 1:nmoths, 'UniformOutput', false));
        lgd.BoxFace.ColorType = 'truecoloralpha';
        lgd.BoxFace.ColorData = uint8(255*[1 1 1 1]');
    end
    % Default all text on this axis to be larger
    set(ax(row,col), 'FontSize', 12)
end
% Set all to same x,y limits
linkaxes(ax)
% X and Y labels
han=axes(gcf,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('r_c (ms)', 'fontsize', 15)
han.XLabel.Position(2) = -0.05;
ylabel('MI (bits / wing stroke)', 'fontsize', 15)
han.YLabel.Position(1) = -0.08;
% Save figure
exportgraphics(gcf,fullfile('figures','KSG_all.pdf'),'ContentType','vector')




%% All KSG precision values in boxplots
cols = {'#CE912D', '#C8C130', '#B96C93', '#0A9167', '#1A689E'...
    '#CE912D', '#C8C130', '#B96C93', '#0A9167', '#1A689E'};
muscle_names = {'L3AX'; 'LBA'; 'LSA'; 'LDVM'; 'LDLM'; 'R3AX'; 'RBA'; 'RSA'; 'RDVM'; 'RDLM'};
jitter_amp = 0.2;
use_fontsize = 15;

%--- All muscles individually
figure('OuterPosition', [844 287 1000 500])
hold on
box on
for i = 1:nmuscles
    boxchart(i*ones(nmoths,1), precision(:,i),...
        'BoxFaceColor', cols{i}, 'BoxFaceAlpha', 0.4,...
        'MarkerColor', 'r', 'WhiskerLineColor', 'k')
    % Get outliers so they don't get jittered
    outlier = isoutlier(precision(:,i), 'quartiles');
    plot(i*ones(sum(~outlier),1) - jitter_amp/2 + jitter_amp*rand(sum(~outlier), 1), precision(~outlier,i), ...
        'ko', 'Markersize', 5, 'MarkerFaceColor', 'w')
    plot(i*ones(sum(outlier),1), precision(outlier,i), ...
        'ko', 'Markersize', 5, 'MarkerFaceColor', 'w')
end
ax = gca;
set(ax, 'XTick', 1:nmuscles)
xt = ax.XTick;
xticklabels([])
for i = 1:length(xt)
    text(xt(i), 0, muscle_names{i}, 'Color', cols{i}, ...
        'Horiz','center', 'Vert','top', 'FontSize', use_fontsize, 'FontWeight', 'bold') 
end
ax.FontSize = use_fontsize;
xlim([0,11])
ylabel('Spike Timing Precision r_c (ms)', 'FontSize', use_fontsize)
% Save figure
exportgraphics(gcf,fullfile('figures','KSG_boxplots_leftright.pdf'),'ContentType','vector')

% Run statistical tests
p_anova = anova1(precision, [], 'off');
[p_kw, tbl, stats] = kruskalwallis(precision, [], 'off');
p_anova
p_kw

%--- Muscles, left+right together
precision_both = [precision(:,1:5); precision(:,6:end)];
figure('OuterPosition', [844 287 450 500])
hold on
box on
for i = 1:5
    boxchart(i*ones(nmoths*2,1), precision_both(:,i),...
        'BoxFaceColor', cols{i}, 'BoxFaceAlpha', 0.4,...
        'MarkerColor', 'r', 'WhiskerLineColor', 'k')
    % Get outliers so they don't get jittered
    outlier = isoutlier(precision_both(:,i), 'quartiles');
    plot(i*ones(sum(~outlier),1) - jitter_amp/2 + jitter_amp*rand(sum(~outlier), 1), precision_both(~outlier,i), ...
        'ko', 'Markersize', 5, 'MarkerFaceColor', 'w')
    plot(i*ones(sum(outlier),1), precision_both(outlier,i), ...
        'ko', 'Markersize', 5, 'MarkerFaceColor', 'w')
end
ax = gca;
set(ax, 'XTick', 1:5)
xt = ax.XTick;
xticklabels([])
for i = 1:length(xt)
    text(xt(i), 0, muscle_names{i}(2:end), 'Color', cols{i}, ...
        'Horiz','center', 'Vert','top', 'Fontsize', use_fontsize, 'FontWeight', 'bold')
end
ax.FontSize = use_fontsize;
xlim([0,6])
ylabel('Spike Timing Precision r_c (ms)')
% Save figure
exportgraphics(gcf,fullfile('figures','KSG_boxplots_combined.pdf'),'ContentType','vector')

% Run statistical tests
p_anova = anova1(precision_both, [], 'off');
[p_kw, tbl, stats] = kruskalwallis(precision_both, [], 'off');
p_anova
p_kw
figure
multcompare(stats);


%--- Spike count
mean_spike_count = zeros(nmoths, nmuscles);
% Preallocate 
for i = 1:nmoths
    load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
    fields = fieldnames(time_data);
    for j = 1:nmuscles
        mean_spike_count(i,j) = mean(sum(~isnan(time_data.(fields{j})), 2));
    end
end
mean_spike_count = [mean_spike_count(:,1:5); mean_spike_count(:,6:end)];
figure('OuterPosition', [844 287 450 500])
hold on
box on
for i = 1:5
    boxchart(i*ones(nmoths*2,1), mean_spike_count(:,i),...
        'BoxFaceColor', cols{i}, 'BoxFaceAlpha', 0.4,...
        'MarkerColor', 'r', 'WhiskerLineColor', 'k');
    % Get outliers so they don't get jittered
    outlier = isoutlier(mean_spike_count(:,i), 'quartiles');
    plot(i*ones(sum(~outlier),1) - jitter_amp/2 + jitter_amp*rand(sum(~outlier), 1), mean_spike_count(~outlier,i), ...
        'ko', 'Markersize', 5, 'MarkerFaceColor', 'w')
    plot(i*ones(sum(outlier),1), mean_spike_count(outlier,i), ...
        'ko', 'Markersize', 5, 'MarkerFaceColor', 'w')
end
ax = gca;
set(ax, 'XTick', 1:5)
xt = ax.XTick;
xticklabels([])
for i = 1:length(xt)
    text(xt(i), 0, muscle_names{i}(2:end), 'Color', cols{i}, ...
        'Horiz','center', 'Vert','top', 'Fontsize', use_fontsize, 'FontWeight', 'bold')
end
ax.FontSize = use_fontsize;
xlim([0,6])
ylabel('Mean Spike Counts per Wing Stroke')
% Save figure
exportgraphics(gcf,fullfile('figures','KSG_boxplots_spike_count.pdf'),'ContentType','vector')

% Run statistical tests
p_anova = anova1(mean_spike_count, [], 'off');
[p_kw, tbl, stats] = kruskalwallis(mean_spike_count, [], 'off');
p_anova
p_kw
figure
multcompare(stats);

%% NSB examples for comparison (Figure 4)
load('NSB_data.mat')
load('NSB_bias_data.mat')

col = '#4472C4';
nsbcol = '#913634';

% Moth 1 L3AX, Moth 3 RSA, Moth 5 LDLM
usemoth = [1, 2, 5];
usemuscle = [1, 8, 5];
for ii = 1:3
    i = usemoth(ii);
    j = usemuscle(ii);
    figure('OuterPosition', [745 496 500 370])
    hold on
    load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
    fields = fieldnames(time_data);
    ind = (i-1)*10 + j;
    % Get bin sizes
    bins = range(time_data.(fields{j}), 'all') ./ (1:nspikebins);
    STD = sqrt(dS_nsbwordvec.^2 + conditionaldS_nsbvec(:, :, ntorque).^2);
    bias = mean(conditionalentropyvec_bias, 4);
    % Plot 
    y_bias = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
    y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque);
    mseb(bins, y, STD(ind, :), struct('col', {{nsbcol}}), 1);
    plot(bins, y_bias, ':', 'color', nsbcol, 'LineWidth', 3, 'DisplayName', 'no')
    % Plot precision dots (in own loop to plot on top of everthing)
    y_bias = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
    compvalue = mean(y_bias(end-10:end));
    compstd = std(y_bias(end-10:end));
    maxval = max(y_bias(y_bias>=compvalue));
    % Peak case, find max value location
    if maxval/compvalue >= 1.2
        [~,precision_ind] = max(y_bias);
    % Plateau case, find farthest right value near compvalue
    else
        precision_ind = find(y_bias >= (compvalue-2*compstd), 1);
    end
    plot(bins(precision_ind), y_bias(precision_ind), '.', 'MarkerSize', 30, 'Color', nsbcol, ...
        'DisplayName', 'no')
    plot(bins(precision_ind), y_bias(precision_ind), 'k.', 'MarkerSize', 15, 'DisplayName', 'no')
    text(0.1, 0.92, ['r_d = ', num2str(round(bins(precision_ind), 3)), ' ms'], ...
            'FontSize', 18, 'units', 'normalized')
    xlim([10^-1.5, 10^2])
    ylim([0, inf])
    set(gca, 'Xscale', 'log')
    xlabel('r_d (ms)', 'Fontsize', 14)
    ylabel('MI (bits / wing stroke)', 'Fontsize', 14)
    set(gca, 'Fontsize', 17)
    exportgraphics(gcf,fullfile('figures',['comparison_NSB_',num2str(ii),'.pdf']),'ContentType','vector')
end

% KSG examples for comparison (figure 4)
usemoth = [1, 2, 5];
usemuscle = [1, 8, 5];
for ii = 1:3
    i = usemoth(ii);
    j = usemuscle(ii);
    figure('OuterPosition', [745 496 500 370])
    hold on
    load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
    fields = fieldnames(time_data);
    meanMI = mean(MI{i,j}, 2) + MI_count(i,j);
    mis = MI_KSG_subsampling_multispike(time_data.(fields{j}), Tz_WSd, knn, (1:nsubsets));
    mi_sd = findMI_KSG_stddev(mis, size(time_data.(fields{j}), 1), false);
    precision_ind = find(meanMI < ((MI{i,j}(1,1) + MI_count(i,j) - mi_sd)), 1);
%     precision_ind = find(meanMI < ((MI{i,j}(1,1) - mi_sd)), 1);
    mseb(noise, meanMI, std(MI{i,j}, 0, 2)', struct('col', {{'#2D427E'}}), 1);
    rectangle('Position', [10^-1.5, meanMI(1)-mi_sd, 10^2-10^-1.5, 2*mi_sd], ...
         'FaceColor', [177 122 168 128] / 256, 'EdgeColor', 'none')
    plot(noise(precision_ind), meanMI(precision_ind), 'r.', 'MarkerSize', 30, 'Color', '#5380FF')
    plot(noise(precision_ind), meanMI(precision_ind), 'k.', 'MarkerSize', 15)
    text(0.16, 0.6, ['r_c = ', num2str(round(noise(precision_ind), 3)), ' ms'], ...
        'FontSize', 18, 'units', 'normalized')
    xlim([10^-1.5, 10^2])
    set(gca, 'XScale', 'log')
    xlabel('r_c (ms)', 'Fontsize', 14)
    ylabel('MI (bits / wing stroke)', 'Fontsize', 14)
    set(gca, 'Fontsize', 17)
    % Save figure
    exportgraphics(gcf,fullfile('figures',['comparison_KSG_',num2str(ii),'.pdf']),'ContentType','vector')
end


%% Compare KSG vs NSB precision values

figure('OuterPosition', [517 420 1300 350])
hold on
col = '#4472C4';
nsbcol = '#913634';
offset = 0.03;

% prec_c = [precision(:,1:5); precision(:,6:end)];
% prec_d = [precision_d(:,1:5); precision_d(:,6:end)];
prec_c = precision;
prec_d = precision_d;

x = linspace(1, 3, 10);
errorbar(x-offset, mean(prec_c, 1), std(prec_c, 1), 'o', ...
        'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
        'MarkerSize', 7, 'LineWidth', 2, 'CapSize', 10)
errorbar(x+offset, mean(prec_d, 1), std(prec_d, 1), 'o', ...
        'color', nsbcol, 'Marker', '^', 'MarkerEdgeColor', nsbcol, 'MarkerFaceColor', nsbcol, ...
        'MarkerSize', 7, 'LineWidth', 2, 'CapSize', 10)
for i = 1:10
    plot(x(i)-offset + 0.1 * normrnd(0, 0.1, 1, nmoths), prec_c(:,i), 'k.')
    plot(x(i)+offset + 0.1 * normrnd(0, 0.1, 1, nmoths), prec_d(:,i), 'k.')
end
xlim([x(1)-0.1, x(end)+0.1])
ylim([0, 13])
set(gca, 'Xtick', x)
% set(gca, 'XTickLabel', cellfun(@(x) x(2:end), muscle_names(1:5), 'UniformOutput', false))
set(gca, 'XTickLabel', muscle_names, 'Fontsize', 14)
ylabel('Precision (ms)', 'Fontsize', 14)
exportgraphics(gcf,fullfile('figures','comparison_KSG_NSB_precision.pdf'),'ContentType','vector')

for i = 1:nmuscles
    ttest2(prec_c(:,i), prec_d(:,i))
end