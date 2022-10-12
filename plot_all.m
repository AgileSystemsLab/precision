% Plot all KSG and NSB

% Load data
load('KSG_data.mat')
load('NSB_data.mat')
load('NSB_bias_data.mat')

nmoths = 7;
nmuscles = 10;
% NSB parameters
ntorque = 2;
nspikebins = 70;

% KSG 2nd derivative peak method
min_peak_height = 0.23535;
min_peak_dist = 10;
sg_ord = 2; % Savitsky-golay filter order
sg_window = 11; % Savitsky-golay filter window length
[b,g] = sgolay(sg_ord, sg_window);

%% Plot all NSB results
STD = sqrt(dS_nsbwordvec.^2 + conditionaldS_nsbvec(:, :, ntorque).^2);
bias = mean(conditionalentropyvec_bias, 4);

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
        x = log10(bins);
        y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
        set(gcf, 'CurrentAxes', ax(row,col))
        mseb(x, y, STD(ind, :), struct('col', {{cols{i}}}), 1);

%         % Get smooth derivative of MI 
%         % (Note magnitude of derivative is currently arbitrary)
%         [b,g] = sgolay(3, sgwindow);
%         ypad = [repmat(y(1),1,sgwindow), y, repmat(y(end),1,sgwindow)];
%         grad = conv(ypad, -1 * g(:,2), 'same'); 
%         grad = grad(sgwindow+1:end-sgwindow+1);
%         plot(grad, 'color', cols{i})

%         plot(x, gradient(y) ./ gradient(x), 'color', cols{i})


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
    % Remove xticks on all but bottom
    if row~=5
        ax(row,col).set('xticklabels', []);
    end
    % Remove yticks on right column
    if col==2
        ax(row,col).set('yticklabels', []);
    end
end
% Set all to same x,y limits
linkaxes(ax)
% X and Y labels
han=axes(gcf,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('log_1_0(r_d (ms))', 'fontsize', 15)
han.XLabel.Position(2) = -0.05;
ylabel('MI (bits / wing stroke)', 'fontsize', 15)
han.YLabel.Position(1) = -0.08;

%% NSB examples for comparison (Figure 4)

% Moth 1 L3AX
figure
hold on
i = 1; j = 1;
ind = (i-1)*10 + j;
load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
fields = fieldnames(time_data);
bins = range(time_data.(fields{j}), 'all') ./ (1:nspikebins);
x = bins;
y_nobias = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque);
y_bias = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
mseb(x, y_nobias, STD(ind, :), struct('col', {{'#2D427E'}}), 1);
mseb(x, y_bias, STD(ind, :), struct('col', {{'#913634'}}), 1);

xlim([0.04, 11])
set(gca, 'XScale', 'log')
set(gca, 'xticklabel', arrayfun(@(x) num2str(x), round(get(gca, 'xtick'), 1), 'UniformOutput', false))


% Moth 5 LDLM
i = 5; j = 5;
ind = (i-1)*10 + j;
figure
hold on
load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
fields = fieldnames(time_data);
bins = range(time_data.(fields{j}), 'all') ./ (1:nspikebins);
x = bins;
y_nobias = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque);
y_bias = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
mseb(x, y_nobias, STD(ind, :), struct('col', {{'#2D427E'}}), 1);
mseb(x, y_bias, STD(ind, :), struct('col', {{'#913634'}}), 1);

xlim([0.04, 11])
set(gca, 'XScale', 'log')
set(gca, 'xticklabel', arrayfun(@(x) num2str(x), round(get(gca, 'xtick'), 1), 'UniformOutput', false))


%% KSG examples for comparison (figure 4)

% Moth 1 L3AX
figure
hold on
i = 1; j = 1;
meanMI = mean(MI{i,j}, 2);

pad = [nan(1, sg_window), meanMI(2)', meanMI(2:end)', nan(1, sg_window)];
grad = conv(pad, -1 * g(:,2), 'same'); 
grad = grad(sg_window+1:end-sg_window);
dpad = [nan(1, sg_window), grad, nan(1, sg_window)];
dgrad = conv(dpad, -1 * g(:,2), 'same'); 
dgrad = dgrad(sg_window+1:end-sg_window);
[~,inds] = findpeaks(dgrad/min(dgrad), 'MinPeakHeight', min_peak_height, 'MinPeakDistance', 20);
if isempty(inds)
    [~,ind] = min(dgrad);
    precision = noise(ind);
else
    ind = inds(1);
    precision = noise(ind);
end

yyaxis right
hold on
xh = xline(precision, 'k--', 'LineWidth', 3);
% Set ConstantLine to 'back'
% Temporarily disable the warning that appears when accessing 
% the undocumented properties.
warnState = warning('off','MATLAB:structOnObject');
cleanupObj = onCleanup(@()warning(warnState)); 
Sxh = struct(xh);        % Get undocumented properties (you'll get a warning)
clear('cleanupObj')      % Trigger warning reset
Sxh.Edge.Layer = 'back'; % Set ConstantLine uistack

plot(noise, -dgrad/min(dgrad))
% dot = plot(precision, -dgrad(ind)/min(dgrad), 'ro', 'MarkerSize', 20, 'MarkerFaceColor', 'r');
ylim([-1.1, 0])

yyaxis left 
mseb(noise, mean(MI{i,j}, 2), std(MI{i,j}, 0, 2)', ...
    struct('col', {{'#2D427E'}}), 1);
xlim([0.04, 11])
set(gca, 'XScale', 'log')
set(gca, 'xticklabel', arrayfun(@(x) num2str(x), round(get(gca, 'xtick'), 1), 'UniformOutput', false))

yyaxis right
radius = 0.02;
xl = get(gca, 'xlim');
yl = range(get(gca, 'ylim'));
x1 = log10(precision);
x2 = x1 + 0.8*radius*range(log10(xl));
w = 10^x2 - 10^x1;
rectangle('Position', [precision-w/2, -1-radius*yl/2, w, radius*yl], 'Curvature', [1 1], 'edgecolor', 'red', 'facecolor', 'red')


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
%         mseb(log10(noise), mean(MI{i,j}, 2), std(MI{i,j}, 0, 2)',...
%             struct('col', {{cols{i}}}), 1);
        plot(log10(noise), gradient(mean(MI{i,j}, 2)) ./ gradient(log10(noise)), 'color', cols{i})
    end
end
% % Plot precision points
% precision = zeros(nmoths, nmuscles);
% ksg_smooth = zeros(nmoths, nmuscles);
% for i = 1:nmoths
%     % Load data (for finding precision values)
%     load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
%     for j = 1:nmuscles
%         % Get row and column, set that ax to gca
%         row = mod(j-1,5)+1;
%         col = int8(j>5)+1;
%         set(gcf, 'CurrentAxes', ax(row,col))
%         % Get precision and plot
%         meanMI = mean(MI{i,j}, 2);
%         mis = MI_KSG_subsampling_multispike(time_data.(fields{j}), Tz_WSd, knn, (1:5));
%         mi_sd = findMI_KSG_stddev(mis, size(Tz_WSd,1), false);
%         precision_ind = find(meanMI < ((MI{i,j}(1,1) - mi_sd)), 1);
%         plot(log10(noise(precision_ind)), meanMI(precision_ind), ...
%             '.', 'color', cols{i}, 'MarkerSize', 18)
%         plot(log10(noise(precision_ind)), meanMI(precision_ind), 'k.', 'MarkerSize', 10)
%         precision(i,j) = noise(precision_ind);
%         r = corrcoef(meanMI(2:end), meanMI(1:end-1));
%         ksg_smooth(i,j) = r(1,2);
%     end
% end
% Plot labels and aesthetics
for j = 1:nmuscles
    % Get row and column, set that ax to gca
    row = mod(j-1,5)+1;
    col = int8(j>5)+1;
    set(gcf, 'CurrentAxes', ax(row,col))
    % Subplot title
    text(0.1, 0.9, fields{j}(1:end-7), 'units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold')
    % Remove xticks on all but bottom
    if row~=5
        ax(row,col).set('xticklabels', []);
    end
    % Remove yticks on right column
    if col==2
        ax(row,col).set('yticklabels', []);
    end
end
% Set all to same x,y limits
linkaxes(ax)
% X and Y labels
han=axes(gcf,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel('log_1_0(r_c (ms))', 'fontsize', 15)
han.XLabel.Position(2) = -0.05;
ylabel('MI (bits / wing stroke)', 'fontsize', 15)
han.YLabel.Position(1) = -0.08;