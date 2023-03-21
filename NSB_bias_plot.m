% NSB Shuffling Bias Correction Supplemental Figure
load('NSB_data.mat')
load('NSB_bias_data.mat')
nmoths = 7;
nmuscles = 10;
% NSB parameters
nspikebins = 70;
% Moth example to use for each muscle (randomly chosen)
usemoth = [6, 3, 2, 4, 1; 7, 5, 1, 3, 6];
muscle_names = {'L3AX'; 'LBA'; 'LSA'; 'LDVM'; 'LDLM'; 'R3AX'; 'RBA'; 'RSA'; 'RDVM'; 'RDLM'};

% Setup main figure
figure('Outerposition', [413, 417, 1017, 456])
ax = gobjects(2, 5);
for i = 1:2
    for j = 1:5
        ax(i,j) = subaxis(2, 5, j, i, 'SpacingVert', 0.03, 'SpacingHoriz', 0.01);
        hold on
    end
end

for j = 1:nmuscles
    col = mod(j-1,5)+1;
    row = int8(j>5)+1;
    i = usemoth(row, col);
    load(fullfile('SubmittedDataallmusclesAllareTzWsd', ['Moth', num2str(i), '_MIdata.mat']))
    fields = fieldnames(time_data);
    ind = (i-1)*10 + j;
    cols = copper(size(conditionaldS_nsbvec, 3));
    set(gcf, 'CurrentAxes', ax(row,col))
    % Get bin sizes
    bins = range(time_data.(fields{j}), 'all') ./ (1:nspikebins);
    % Plot MI vs precision at different ntorques
    for ntorque = 2:size(conditionaldS_nsbvec, 3)
        STD = sqrt(dS_nsbwordvec.^2 + conditionaldS_nsbvec(:, :, ntorque).^2);
        bias = mean(conditionalentropyvec_bias, 4);
        % Plot 
        y_bias = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque) - (S_nsbwordvec(ind,:) - bias(ind,:,ntorque));
        y = S_nsbwordvec(ind, :) - conditionalentropyvec(ind, :, ntorque);
        mseb(bins, y, STD(ind, :), struct('col', {{cols(ntorque,:)}}), 1);
        plot(bins, y_bias, ':', 'color', cols(ntorque,:), 'LineWidth', 3, 'DisplayName', 'no')
    end
    % Plot precision dots (in own loop to plot on top of everthing)
    for ntorque = 2:size(conditionaldS_nsbvec, 3)
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
        plot(bins(precision_ind), y_bias(precision_ind), '.', 'MarkerSize', 30, 'Color', cols(ntorque,:), ...
            'DisplayName', 'no')
        plot(bins(precision_ind), y_bias(precision_ind), 'k.', 'MarkerSize', 15, 'DisplayName', 'no')
    end
    text(0.02, 0.93, [muscle_names{j}, ', Moth ', num2str(i)], ...
        'FontSize', 14, 'units', 'normalized')
    xlim([10^-1.5, 10^2])
    ylim([0, 3])
    set(gca, 'Xscale', 'log')
    if row == 2
        xlabel('r_d (ms)')
    else
        set(gca, 'XTickLabel', [])
    end
    if col == 1
        ylabel('MI (bits / wing stroke)')
    else
        set(gca, 'YTickLabel', [])
    end
    if col == 5
        set_leg_off = findobj('DisplayName', 'no');
        for k = 1:numel(set_leg_off)
            set_leg_off(k).Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        lgd = legend({'4', '9', '16', '25'});
        title(lgd, '# of motor states')
        lgd.BoxFace.ColorType='truecoloralpha';
        lgd.BoxFace.ColorData=uint8(255*[1 1 1 0.25]');
        lgd.Position = lgd.Position - [0 0.05 0 0];
    end
end

% Save figure
exportgraphics(gcf, fullfile('figures','NSB_bias_supp.pdf'),'ContentType','vector')