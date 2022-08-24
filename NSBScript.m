% 

nmoths = 7;
nspikebins = 70;
ntorque = 3;
sgwindow = 13;
deriv_thresh = 0.04;

% Load NSB results
load('NSB_data.mat')
load('NSB_bias_data.mat')


% Plot all NSB results
STD = sqrt(dS_nsbwordvec.^2 + conditionaldS_nsbvec(:, :, ntorque).^2);
bias = mean(conditionalentropyvec_bias, 4);

% Set up figure
figure('Outerposition', [597, 61, 833, 812])
ax = gobjects(5, 2);
for i = 1:5
    for j = 1:2
        ax(i,j) = subaxis(5, 2, j, i, 'SpacingVert', 0.01);
        hold on
    end
end
cols = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
precision = nan(nmoths, nmuscles);
precision_y = nan(nmoths, nmuscles);
smoothness = zeros(nmoths, nmuscles);
% Loop over moths
for i = 1:nmoths
    load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
    fields = fieldnames(time_data);
    % Loop over muscles
    for j = 1:length(fields)
        % Get bin sizes
        bins = (max(time_data.(fields{j}), [], 'all') - min(time_data.(fields{j}), [], 'all')) ./ (1:nspikebins);
        % Plot 
        row = mod(j-1,5)+1;
        col = int8(j>5)+1;
        ind = (i-1)*10 + j;
        % Plot main precision line
        x = bins;
%         y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque);
        y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
%         y = S_nsbwordvec(ind,:) - bias(ind,:,ntorque);
        set(gcf, 'CurrentAxes', ax(row,col))
        mseb(x, y, STD(ind, :), struct('col', {{cols{i}}}), 1);
        if i==1 && j==1
            title([num2str(ntorque),' torque bins'])
        end

        % Get smooth derivative of MI 
        % (Note units of derivative are arbitrary)
        [b,g] = sgolay(3, sgwindow);
        ypad = [repmat(y(1),1,sgwindow), y, repmat(y(end),1,sgwindow)];
        grad = conv(ypad, -1 * g(:,2), 'same'); 
        grad = grad(sgwindow+1:end-sgwindow+1);
%         plot(grad, 'color', cols{i})

        % Precision is first zero-crossing of derivative
        % IF derivative gets above threshold before passing zero
        flip = find(sign(grad)==-1, 1);
        if (flip~=1) && all(grad(1:flip-1)>0) && (max(grad(1:flip-1)) >= deriv_thresh)
            precision_ind = flip-1;
            precision(i,j) = log10(bins(precision_ind));
            precision_y(i,j) = y(precision_ind);
        end
        
        % Quantify smoothness metric
%         smoothness(i,j) = trapz(diff(diff(y)).^2);
        r = corrcoef(y(2:end), y(1:end-1));
        smoothness(i,j) = r(1,2);

%         % Segmented regression to find precision
%         lerr = zeros(nspikebins,1); 
%         rerr = zeros(nspikebins,1);
%         x = (1:nspikebins);
%         % Loop over segment points (s is shared)
%         for s = 3:(size(S_nsbwordvec, 2) - 2)
%             % Left orthogonal regression line fit
%             v = pca([x(1:s)' y(1:s)']);
%             beta = v(2,1)/v(1,1);
%             int = mean(y(1:s)) - beta * mean(x(1:s));
%             lerr(s) = pointLineSSerr(x(1:s), y(1:s), [beta, int]);
%             % Right orthogonal regression line fit
%             v = pca([x(s:end)' y(s:end)']);
%             beta = v(2,1)/v(1,1);
%             int = mean(y(s:end)) - beta * mean(x(s:end));
%             rerr(s) = pointLineSSerr(x(s:end), y(s:end), [beta, int]);
%         end
%         err = lerr + rerr;
%         [val,precision] = min(err(err~=0));
%         minerr(i,j) = val;
%         precision = precision + 2;
%         plot(ax(row,col), log10(bins(precision)), y(precision), 'r*')

    end
end
for i = 1:nmoths
    for j = 1:nmuscles
        row = mod(j-1,5)+1;
        col = int8(j>5)+1;
        set(gcf, 'CurrentAxes', ax(row,col))
        plot(precision(i,j), precision_y(i,j), ...
                '.', 'color', cols{i}, 'MarkerSize', 18)
        plot(precision(i,j), precision_y(i,j), 'k.', 'MarkerSize', 10)
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

figure('Outerposition', [597, 61, 833, 812])
ax = gobjects(5, 2);
for i = 1:5
    for j = 1:2
        ax(i,j) = subaxis(5, 2, j, i, 'SpacingVert', 0.01);
        hold on
    end
end
% Loop over moths
for i = 1:nmoths
    load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
    fields = fieldnames(time_data);
    % Loop over muscles
    for j = 1:length(fields)
        % Get bin sizes
        bins = (max(time_data.(fields{j}), [], 'all') - min(time_data.(fields{j}), [], 'all')) ./ (1:nspikebins);
        % Plot 
        row = mod(j-1,5)+1;
        col = int8(j>5)+1;
        ind = (i-1)*10 + j;
        % Plot main precision line
%         y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque);
        y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
        set(gcf, 'CurrentAxes', ax(row,col))
        % Get smooth derivative of MI 
        % (Note units of derivative are arbitrary)
        [b,g] = sgolay(3, sgwindow);
        ypad = [repmat(y(1),1,sgwindow), y, repmat(y(end),1,sgwindow)];
        grad = conv(ypad, -1 * g(:,2), 'same'); 
        grad = grad(sgwindow+1:end-sgwindow);
        plot(grad, 'color', cols{i})
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


%% Precision boxplots
figure
boxplot(10 .^ precision, 'Labels', cellfun(@(x) x(1:end-7), fields, 'UniformOutput', false))
title([num2str(ntorque),' torque bins'])

% %% Segmented Regression Test
% 
% i = 5; % moth
% j = 5; % muscle
% ind = (i-1)*10 + j;
% x = (1:nspikebins);
% y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, 2);
% 
% % Preallocate
% lbeta = zeros(nspikebins, 2);
% rbeta = zeros(nspikebins, 2);
% lerr = zeros(nspikebins, 1);
% rerr = zeros(nspikebins, 1);
% % Loop over segment points (s is shared)
% for s = 3:(size(S_nsbwordvec, 2) - 2)
%     % Left orthogonal regression line fit
%     v = pca([x(1:s)' y(1:s)']);
%     beta = v(2,1)/v(1,1);
%     int = mean(y(1:s)) - beta * mean(x(1:s));
%     lbeta(s,:) = [beta, int];
%     lerr(s) = pointLineSSerr(x(1:s), y(1:s), lbeta(s,:));
%     % Right orthogonal regression line fit
%     v = pca([x(s:end)' y(s:end)']);
%     beta = v(2,1)/v(1,1);
%     int = mean(y(s:end)) - beta * mean(x(s:end));
%     rbeta(s,:) = [beta, int];
%     rerr(s) = pointLineSSerr(x(s:end), y(s:end), rbeta(s,:));
% end
% 
% figure
% 
% subplot(2,1,1)
% hold on
% plot(lerr)
% plot(rerr)
% plot(lerr+rerr)
% bob = lerr+rerr;
% [val,ind] = min(bob(bob~=0));
% plot(ind+2, val, 'r*')
% title(['Moth ', num2str(i), ' ', fields{j}(1:end-7)])
% 
% subplot(2,1,2)
% hold on
% plot(x,y)
% plot(x(ind+2), y(ind+2), 'r*')

