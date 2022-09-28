% KSG precision results comparison for different precision-picking methods
load('KSG_data.mat')

nmoths = 7;
nmuscles = 10;

savefigs = true;

deriv_thresh_scale = 0.38687; % Derivative threshold
sg_ord = 2; % Savitsky-golay filter order
sg_window = 11; % Savitsky-golay filter window length
[b,g] = sgolay(sg_ord, sg_window);

s = 30; % how many samples on each side to fit line to


% Example of normal method vs. derivative method
examp_moth = 4;
examp_muscle = 4;

load(fullfile('Data',['Moth',num2str(examp_moth),'_MIdata.mat']))
fields = fieldnames(time_data);
mis = MI_KSG_subsampling_multispike(time_data.(fields{examp_muscle}), Tz_WSd, knn, (1:10));
mi_sd = findMI_KSG_stddev(mis, size(time_data.(fields{examp_muscle}), 1), false);
meanMI = mean(MI{examp_moth, examp_muscle}, 2);

figure('outerposition', [440 366 967 350])
t = tiledlayout(2,3);
%--- Original method
precision_ind = find(meanMI < ((meanMI(1) - mi_sd)), 1);
nexttile([2,1])
hold on
rectangle('Position', [log10(noise(2)), meanMI(1) - mi_sd, log10(noise(end))-log10(noise(2)), 2*mi_sd], ...
    'FaceColor', [177 122 168]./256, 'EdgeColor', 'none')
mseb(log10(noise), meanMI, std(MI{examp_moth,examp_muscle}, 0, 2)', struct('col', {{'#2D427E'}}), 1);
plot(log10(noise(precision_ind)), meanMI(precision_ind), 'r.', 'MarkerSize', 15)
ylim([0, 0.9])
xlabel('log_1_0(r_c (ms))')
ylabel('MI (bits / wing stroke)')
title('STD Threshold Method')

%--- Derivative method
pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
grad = conv(pad, -1 * g(:,2), 'same'); 
grad = grad(sg_window+1:end-sg_window);
deriv_thresh = min(grad) * deriv_thresh_scale;
deriv_prec_ind = find(grad < deriv_thresh, 1);
% regular MI vs noise
nexttile([1,1])
hold on
mseb(log10(noise), meanMI, std(MI{examp_moth,examp_muscle}, 0, 2)', struct('col', {{'#2D427E'}}), 1);
plot(log10(noise(deriv_prec_ind)), meanMI(deriv_prec_ind), 'r.', 'MarkerSize', 15)
ylim([-inf 0.9])
xlabel('log_1_0(r_c (ms))')
ylabel('MI (bits/ws)')
title('Derivative Method')


%--- Two-line Method
x = log10(noise);
v = pca([x(2:s+1)' meanMI(2:s+1)]);
beta = v(2,1)/v(1,1);
int = mean(meanMI(2:s+1)) - beta * mean(x(2:s+1));
lbeta = [beta, int];
v = pca([x(end-s:end)' meanMI(end-s:end)]);
beta = v(2,1)/v(1,1);
int = mean(meanMI(end-s:end)) - beta * mean(x(end-s:end));
rbeta = [beta, int];
twoline_precision = (rbeta(2) - lbeta(2)) / (lbeta(1) - rbeta(1));
% Plot 
nexttile([2,1])
hold on
plot(log10(noise), lbeta(1) * log10(noise) + lbeta(2), 'k--', 'LineWidth', 2)
plot(log10(noise), rbeta(1) * log10(noise) + rbeta(2), 'k--', 'LineWidth', 2)
mseb(log10(noise), meanMI, std(MI{examp_moth,examp_muscle}, 0, 2)', struct('col', {{'#2D427E'}}), 1);
plot(log10(noise(2:s+1)), meanMI(2:s+1), '.-', 'Color', '#F86E52')
plot(log10(noise(end-s:end)), meanMI(end-s:end), '.-', 'Color', '#F86E52')
plot(twoline_precision, rbeta(1) * twoline_precision + rbeta(2), 'r.', 'MarkerSize', 15)
ylim([0, 0.9])
xlabel('log_1_0(r_c (ms))')
ylabel('MI (bits/ws)')
title('Line Intersection Method')

% Derivative of MI vs noise (down here because nexttile is dumb)
nexttile([1,1])
hold on
plot(log10(noise), -grad/min(grad), 'Color', '#2D427E')
yline(-deriv_thresh_scale, 'LineWidth', 2)
plot(log10(noise(deriv_prec_ind)), -grad(deriv_prec_ind)/min(grad), 'r.', 'MarkerSize', 15)
xlabel('log_1_0(r_c (ms))')
ylabel('d MI / d r_c (a.u.)')

if savefigs
    exportgraphics(gcf,fullfile('figures','KSG_precision_methods.pdf'),'ContentType','vector')
end

% Precision values
% Original method
precision = zeros(nmoths, nmuscles);
deriv_precision = zeros(nmoths, nmuscles);
twoline_precision = zeros(nmoths, nmuscles);
x = log10(noise);
for i = 1:nmoths
    % Load data
    load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
    fields = fieldnames(time_data);
    for j = 1:nmuscles
        meanMI = mean(MI{i,j}, 2);
        mis = MI_KSG_subsampling_multispike(time_data.(fields{j}), Tz_WSd, knn, (1:10));
        mi_sd = findMI_KSG_stddev(mis, size(time_data.(fields{j}), 1), false);
        precision_ind = find(meanMI < ((MI{i,j}(1,1) - mi_sd)), 1);
        precision(i,j) = noise(precision_ind);
        % Alternative ways to find precision
        % Get smooth derivative of MI, find where passes threshold
        pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
        grad = conv(pad, -1 * g(:,2), 'same'); 
        grad = grad(sg_window+1:end-sg_window);
        deriv_thresh = min(grad) * deriv_thresh_scale;
        deriv_prec_ind = find(grad < deriv_thresh, 1);
        deriv_precision(i,j) = noise(deriv_prec_ind);
        % Two-line orthogonal regression meeting point
        v = pca([x(2:s+1)' meanMI(2:s+1)]);
        beta = v(2,1)/v(1,1);
        int = mean(meanMI(2:s+1)) - beta * mean(x(2:s+1));
        lbeta = [beta, int];
        v = pca([x(end-s:end)' meanMI(end-s:end)]);
        beta = v(2,1)/v(1,1);
        int = mean(meanMI(end-s:end)) - beta * mean(x(end-s:end));
        rbeta = [beta, int];
        twoline_precision(i,j) = 10^((rbeta(2) - lbeta(2)) / (lbeta(1) - rbeta(1)));
    end
end

figure('Outerposition', [440, 495, 967, 250])
ax = gobjects(1,5);
col = '#4472C4';
h_deriv = zeros(1,5); p_deriv = zeros(1,5);
h_intersect = zeros(1,5); p_intersect = zeros(1,5);
for i = 1:5
    ax(i) = subplot(1,5,i);
    hold on
    X = 0.1 * normrnd(0, 0.5, 1, 2*nmoths);
    dX = 1 + 0.1 * normrnd(0, 0.5, 1, 2*nmoths);
    tX = 2 + 0.1 * normrnd(0, 0.5, 1, 2*nmoths);
    P = [precision(:,i); precision(:,5+i)];
    dP = [deriv_precision(:,i); deriv_precision(:,5+i)];
    tP = [twoline_precision(:,i); twoline_precision(:,5+i)];
    errorbar([0 1 2], [mean(P) mean(dP) mean(tP)], [std(P) std(dP) std(tP)], 'o', ...
        'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
        'LineWidth', 1, 'CapSize', 13)
    plot(X, P, 'k.', 'MarkerSize', 7, 'color', [0 0 0 0.3])
    plot(dX, dP, 'r.', 'MarkerSize', 7, 'color', [0 0 0 0.3])
    plot(tX, tP, 'r.', 'MarkerSize', 7, 'color', [0 0 0 0.3])
    
    xlim([-0.5 2.5])
    ylim([0, 4])
    set(gca, 'xtick', [0 1 2]);
    set(gca, 'xticklabels', {'STD', 'Derivative', 'Intersection'});
    title(fields{i}(2:end-7))
    [h_deriv(i), p_deriv(i)] = ttest2(P, dP);
    [h_intersect(i), p_intersect(i)] = ttest2(P, tP);
end
subplot(1,5,1)
ylabel('Precision (ms)')
linkaxes(ax, 'y')

h_deriv
p_deriv
h_intersect
p_intersect
if savefigs
    exportgraphics(gcf,fullfile('figures','KSG_precision_methods_distributions.pdf'),'ContentType','vector')
end