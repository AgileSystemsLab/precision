rng('shuffle')
% Controls
do_long_runs = false;
do_threshold_finding = false;
do_method_comparison = true;

% Main constants
nmoths = 7;
nmuscles = 10;
knn = 4;
noise = [0, logspace(log10(0.05), log10(6), 120)];
repeats = 150;
prec_levels = 1:0.5:4;
n = length(prec_levels);
% Derivative method 
deriv_thresh_scale = 0.38687;
sg_ord = 2; % Savitsky-golay filter order
sg_window = 11; % Savitsky-golay filter window length
[b,g] = sgolay(sg_ord, sg_window);
% Two-line method
s = 5; % how many samples on each side to fit line to
x = log10(noise);


%% Fixed precision on real dataset
if do_long_runs
    % Preallocate
    MI_disc = cell(nmoths, nmuscles, n);
    MI_disc(:) = {zeros(length(noise), repeats)};
    fakeX = cell(nmoths, nmuscles, n);
    % Loop over moths
    for i = 1:nmoths
        % Load data
        load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
        fields = fieldnames(time_data);
        disp(['Moth ',num2str(i)])
        % Loop over muscles
        for j = 1:length(fields)
            disp(fields{j}(1:end-7))
            % Create fixed precision real spike data as "fake" X, run precision estimation
            % Loop over precision levels
            for k = 1:n
                fakeX{i, j, k} = round(time_data.(fields{j}) / prec_levels(k)) * prec_levels(k);
                MI_disc{i, j, k} = KSG_precision(fakeX{i,j,k}, Tz_WSd, knn, repeats, noise);
            end
        end
    end
    % Save results
    save('KSG_sim_real_data_fixed_precision.mat', 'MI_disc', 'fakeX', 'noise', 'repeats', 'knn')
end

%% Fixed precision on fully synthetic dataset from bivariate gaussian
repeats_at_corr = 1;
num_points = 2500;
corr = [0.5, 0.6, 0.7, 0.8, 0.9];
ncorr = length(corr);
if do_long_runs
    muX = 0;         %X mean of the bivariate gaussian
    muY = 0;         %Y mean of the bivariate gaussian
    sigmaX = 2;      %X standard deviation   
    sigmaY = 2;      %Y standard deviation
    
    % preallocate
    MI_synth = cell(ncorr, repeats_at_corr, n);
    MI_synth(:) = {zeros(length(noise), repeats)};
    for i = 1:ncorr
        for r = 1:repeats_at_corr
            %generate a correlated X and Y
            x1 = normrnd(0, 1, num_points, 1);
            x2 = normrnd(0, 1, num_points, 1);
            x3 = corr(i) .* x1 + (1 - corr(i)^2)^.5 .* x2;
            synthX = muX + x1 * sigmaX;
            synthY = muY + x3 * sigmaY;
            % Loop over fixed precision levels
            for k = 1:n
                synthX_disc = round(synthX / prec_levels(k)) * prec_levels(k);
                MI_synth{i,j,k} = KSG_precision(synthX_disc, synthY, knn, repeats, noise);
            end
        end
    end
    
    % Save results
    save('KSG_sim_synthetic_data_fixed_precision.mat', 'MI_synth', ...
        'noise', 'repeats', 'knn', 'corr', 'repeats_at_corr', 'num_points', 'correlation')
end

%% Find optimal derivative threshold
if do_threshold_finding
    load('KSG_sim_real_data_fixed_precision.mat')
    load('KSG_sim_synthetic_data_fixed_precision.mat')
    thresh_scale = linspace(0.2, 0.7, 100);
    %------ Real dataset
    err = zeros(size(thresh_scale));
    % Loop over different derivative thresholds
    for ii = 1:length(thresh_scale)
        deriv_thresh_scale = thresh_scale(ii);
        deriv_precision = cell(nmoths, nmuscles);
        deriv_precision(:) = {nan(1, n)};
        for i = 1:nmoths
            for j = 1:nmuscles
                for k = 1:n
                    pad = [nan(1, sg_window), mean(MI_disc{i,j,k}, 2)', nan(1, sg_window)];
                    grad = conv(pad, -1 * g(:,2), 'same'); 
                    grad = grad(sg_window+1:end-sg_window);
                    deriv_thresh = min(grad) * deriv_thresh_scale;
                    deriv_prec_ind = find(grad < deriv_thresh, 1);
                    deriv_precision{i,j}(k) = noise(deriv_prec_ind);
                end
            end
        end
        err(ii) = sqrt(sum(cellfun(@(x) sum((x - prec_levels).^2), deriv_precision), 'all'));
    end
    % Real data plot 
    figure
    hold on
    plot(thresh_scale, err)
    [~,ind] = min(err);
    plot(thresh_scale(ind), err(ind), 'r.', markersize=25)
    text(0.35, 0.4, ['Optimal threshold = ', num2str(thresh_scale(ind))], 'units', 'normalized')
    title('Real dataset with controlled precision')
    xlabel('Derivative threshold')
    ylabel('Sum of squared errors')
    exportgraphics(gcf,fullfile('figures','threshold_finding_realdataset.pdf'),'ContentType','vector')

    %------ Synthetic dataset
    err = zeros(size(thresh_scale));
    % Loop over different derivative thresholds
    for ii = 1:length(thresh_scale)
        deriv_thresh_scale = thresh_scale(ii);
        synth_precision = zeros(ncorr, repeats_at_corr, n);
        for i = 1:ncorr
            for j = 1:repeats_at_corr
                for k = 1:n
                    pad = [nan(1, sg_window), mean(MI_synth{i,j,k}, 2)', nan(1, sg_window)];
                    grad = conv(pad, -1 * g(:,2), 'same'); 
                    grad = grad(sg_window+1:end-sg_window);
                    deriv_thresh = min(grad) * deriv_thresh_scale;
                    ind = find(grad < deriv_thresh, 1);
                    synth_precision(i,j,k) = noise(ind);
                end
            end
        end
        synth_precision = reshape(synth_precision, [ncorr*repeats_at_corr, n]);
    
        err(ii) = sqrt(sum((synth_precision-prec_levels).^2, 'all'));
    end
    % Synthetic data plot
    figure
    hold on
    plot(thresh_scale, err)
    [~,ind] = min(err);
    plot(thresh_scale(ind), err(ind), 'r.', markersize=25)
    text(0.35, 0.4, ['Optimal threshold = ', num2str(thresh_scale(ind))], 'units', 'normalized')
    title('Synthetic dataset with controlled precision')
    xlabel('Derivative threshold')
    ylabel('Sum of squared errors')
    exportgraphics(gcf,fullfile('figures','threshold_finding_synthdataset.pdf'),'ContentType','vector')
end

%% Methods comparison 
% (really a companion to KSG_precision_comparison but the data for this is
% generated here so kept in this script for simplicity of running)
deriv_thresh_scale = 0.38687;

if do_method_comparison
    load('KSG_sim_synthetic_data_fixed_precision.mat')

    muX = 0;         %X mean of the bivariate gaussian
    muY = 0;         %Y mean of the bivariate gaussian
    sigmaX = 2;      %X standard deviation   
    sigmaY = 2;      %Y standard deviation

    precision_std = zeros(ncorr, repeats_at_corr, n);
    precision_deriv = zeros(ncorr, repeats_at_corr, n);
    precision_twoline = zeros(ncorr, repeats_at_corr, n);
    for i = 1:ncorr
        for j = 1:repeats_at_corr
            %generate a correlated X and Y
            x1 = normrnd(0, 1, num_points, 1);
            x2 = normrnd(0, 1, num_points, 1);
            x3 = corr(i) .* x1 + (1 - corr(i)^2)^.5 .* x2;
            synthX = muX + x1 * sigmaX;
            synthY = muY + x3 * sigmaY;
            for k = 1:n
                meanMI = mean(MI_synth{i,j,k}, 2);
                % STD method
                mis = MI_KSG_subsampling_multispike(synthX, synthY, knn, (1:4));
                mi_sd = findMI_KSG_stddev(mis, size(synthY,1), false);
                ind = find(meanMI < (MI_synth{i,j,k}(1,1) - mi_sd), 1);
                if isempty(ind)
                    precision_std(i,j,k) = nan;
                else
                    precision_std(i,j,k) = noise(ind);
                end
                % Deriv method
                pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
                grad = conv(pad, -1 * g(:,2), 'same'); 
                grad = grad(sg_window+1:end-sg_window);
                deriv_thresh = min(grad) * deriv_thresh_scale;
                ind = find(grad < deriv_thresh, 1);
                precision_deriv(i,j,k) = noise(ind);
                % Two line method
                v = pca([x(2:s+1)' meanMI(2:s+1)]);
                beta = v(2,1)/v(1,1);
                int = mean(meanMI(2:s+1)) - beta * mean(x(2:s+1));
                lbeta = [beta, int];
                v = pca([x(end-s:end)' meanMI(end-s:end)]);
                beta = v(2,1)/v(1,1);
                int = mean(meanMI(end-s:end)) - beta * mean(x(end-s:end));
                rbeta = [beta, int];
                precision_twoline(i,j,k) = 10^((rbeta(2) - lbeta(2)) / (lbeta(1) - rbeta(1)));
            end
        end
    end
    % Reshape for plotting
    precision_std = reshape(precision_std, [ncorr*repeats_at_corr, n]);
    precision_deriv = reshape(precision_deriv, [ncorr*repeats_at_corr, n]);
    precision_twoline = reshape(precision_twoline, [ncorr*repeats_at_corr, n]);
    % Plot
    figure('outerposition', [440 366 967 300])
    t = tiledlayout(1, 3);
    col = '#4472C4';
    % STD method
    nexttile()
    hold on
    errorbar(prec_levels, mean(precision_std, 1), std(precision_std, 1), 'o', ...
        'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
        'LineWidth', 1, 'CapSize', 13)
    xlim([prec_levels(1)-0.25, prec_levels(end)+0.25])
    ylim([0, 5])
    plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
    title('STD Threshold Method')
    xlabel('Actual precision (ms)')
    ylabel('Measured precision (ms)')
    % Derivative method
    nexttile()
    hold on
    errorbar(prec_levels, mean(precision_deriv, 1), std(precision_deriv, 1), 'o', ...
        'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
        'LineWidth', 1, 'CapSize', 13)
    xlim([prec_levels(1)-0.25, prec_levels(end)+0.25])
    ylim([0, 5])
    plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
    title('Derivative Method')
    xlabel('Actual precision (ms)')
    ylabel('Measured precision (ms)')
    % Two line method
    nexttile()
    hold on
    errorbar(prec_levels, mean(precision_twoline, 1), std(precision_twoline, 1), 'o', ...
        'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
        'LineWidth', 1, 'CapSize', 13)
    xlim([prec_levels(1)-0.25, prec_levels(end)+0.25])
    ylim([0, 5])
    plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
    title('Line Intersection Method')
    xlabel('Actual precision (ms)')
    ylabel('Measured precision (ms)')
    % Save 
    exportgraphics(gcf,fullfile('figures','KSG_precision_methods_simulations.pdf'),'ContentType','vector')
end


%% Distributions of precision at specific known levels
% Actually find precision for both real and synthetic fixed datasets
load('KSG_sim_real_data_fixed_precision.mat')
load('KSG_sim_synthetic_data_fixed_precision.mat')
% deriv_thresh_scale = 0.38687;
deriv_thresh_scale = 0.2;
% deriv_thresh = -0.02;
sg_ord = 2; % Savitsky-golay filter order
sg_window = 13; % Savitsky-golay filter window length
[b,g] = sgolay(sg_ord, sg_window);
% Real
figure
hold on
precision_real = zeros(nmoths, nmuscles, n);
for i = 1:nmoths
    for j = 1:nmuscles
        for k = 1:n
            meanMI = mean(MI_disc{i,j,k}, 2);
            % initial MI drop from 0 noise to some noise throws off derivative estimation, smooth out
            meanMI(1) = meanMI(2); 

            pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
            grad = conv(pad, -1 * g(:,2), 'same'); 
            grad = grad(sg_window+1:end-sg_window);

            dpad = [nan(1, sg_window), grad, nan(1, sg_window)];
            dgrad = conv(dpad, -1 * g(:,2), 'same'); 
            dgrad = dgrad(sg_window+1:end-sg_window);
            [~,inds] = findpeaks(dgrad/min(dgrad), 'MinPeakHeight', 0.5);
            if isempty(inds)
                [~,ind] = min(dgrad);
                precision_real(i,j,k) = noise(ind);
            else
                precision_real(i,j,k) = noise(inds(1));
            end
        end
    end
end
% Synthetic
precision_synth = zeros(ncorr, repeats_at_corr, n);
for i = 1:ncorr
    for j = 1:repeats_at_corr
        for k = 1:n
            meanMI = mean(MI_synth{i,j,k}, 2);
            meanMI(1) = meanMI(2); 
            meanMI = meanMI / max(meanMI);
            pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
            grad = conv(pad, -1 * g(:,2), 'same'); 
            grad = grad(sg_window+1:end-sg_window);
            dpad = [nan(1, sg_window), grad, nan(1, sg_window)];
            dgrad = conv(dpad, -1 * g(:,2), 'same'); 
            dgrad = dgrad(sg_window+1:end-sg_window);
%             deriv_thresh = min(grad) * deriv_thresh_scale;
%             ind = find(grad < deriv_thresh, 1);
            [~,inds] = findpeaks(dgrad/min(dgrad), 'MinPeakHeight', 0.3);
            if isempty(inds)
                [~,ind] = min(dgrad);
                precision_synth(i,j,k) = noise(ind);
            else
                precision_synth(i,j,k) = noise(inds(1));
            end
            plot(dgrad/min(dgrad), 'color', cols(k,:))
            if (i == 1) && (j == 1)
                [~,ind] = min(abs(prec_levels(k) - noise));
                xline(ind, 'color', cols(k,:))
            end
        end
    end
end
% Reshape both for ease of plotting
precision_real = reshape(precision_real, [nmoths*nmuscles, n]);
precision_synth = reshape(precision_synth, [ncorr*repeats_at_corr, n]);


figure('OuterPosition', [1007, 234, 280, 610])
col = '#4472C4';
t = tiledlayout(2, 1);

% Real dataset plot
nexttile()
hold on
errorbar(prec_levels, mean(precision_real, 1), std(precision_real, 1), 'o', ...
    'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
    'LineWidth', 1, 'CapSize', 13)
xlim([prec_levels(1)-0.25, prec_levels(end)+0.25])
ylim([0, 5])
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
title('KSG, Real dataset')
xlabel('Actual precision (ms)')
ylabel('Measured precision (ms)')
% Synthetic dataset plot
nexttile()
hold on
errorbar(prec_levels, mean(precision_synth, 1), std(precision_synth, 1), 'o', ...
    'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
    'LineWidth', 1, 'CapSize', 13)
xlim([prec_levels(1)-0.25, prec_levels(end)+0.25])
ylim([0, 5])
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
title('KSG, Synthetic dataset')
xlabel('Actual precision (ms)')
ylabel('Measured precision (ms)')

% exportgraphics(gcf,fullfile('figures','simulations_KSG.pdf'),'ContentType','vector')

%% Example of MI vs noise at different 

load('KSG_data.mat')
deriv_thresh_scale = 0.38687;

i = 3; % moth
j = 3; % muscle

figure 
hold on
cols = copper(n);
mseb(log10(noise), mean(MI{i,j}, 2), std(MI{i,j}, 0, 2)');
for k = 1:n
    meanMI = mean(MI_disc{i,j,k}, 2);
    mseb(log10(noise), meanMI, std(MI_disc{i,j,k}, 0, 2)', ...
        struct('col', {{cols(k,:)}}), 1)
    % Get precision (I know it's stored somewhere but this is faster to write)
    pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
    grad = conv(pad, -1 * g(:,2), 'same'); 
    grad = grad(sg_window+1:end-sg_window);
    deriv_thresh = min(grad) * deriv_thresh_scale;
    ind = find(grad < deriv_thresh, 1);
    plot(log10(noise(ind)), meanMI(ind), '.', 'color', cols(k,:), 'MarkerSize', 25)
    plot(log10(noise(ind)), meanMI(ind), 'k.', 'MarkerSize', 15)
    if k == 1
        bob = grad;
        jim = meanMI;
    end

end


%% Example of synthetic data plot

muX = 0;         %X mean of the bivariate gaussian
muY = 0;         %Y mean of the bivariate gaussian
sigmaX = 2;      %X standard deviation   
sigmaY = 2;      %Y standard deviation
%generate a correlated X and Y
x1 = normrnd(0, 1, num_points, 1);
x2 = normrnd(0, 1, num_points, 1);

figure('outerposition', [841, 455, 700, 350])
ax = gobjects(1, 2);
cols = parula(ncorr+1);
hspacing = 0.15;

ax(1) = subaxis(1, 2 , 1, 'SpacingHoriz', hspacing);
ax(1).Box = 'off';
hold on
for i = 1:ncorr
    x3 = corr(i) .* x1 + (1 - corr(i)^2)^.5 .* x2;
    synthX = muX + x1 * sigmaX;
    synthY = muY + x3 * sigmaY;

    plot(synthX, synthY, '.', 'color', cols(i,:))
end
xlabel('X', 'FontWeight', 'bold', 'FontSize', 12)
ylabel('Y', 'FontWeight', 'bold', 'FontSize', 12)
set(get(gca, 'Ylabel'), 'Rotation', 0)
xlim([-7, 7])
[~,leg] = legend(arrayfun(@(x) ['\rho = ',num2str(x)], corr, 'UniformOutput', false), 'location', 'southeast');
set(findobj(leg, '-property', 'MarkerSize'), 'MarkerSize', 13)

ax(2) = subaxis(1, 2 , 2, 'SpacingHoriz', hspacing);
ax(2).Box = 'off';
hold on
for i = 1:ncorr
    x3 = corr(i) .* x1 + (1 - corr(i)^2)^.5 .* x2;
    synthX = muX + x1 * sigmaX;
    synthY = muY + x3 * sigmaY;

    synthX_disc = round(synthX / 1) * 1;
    plot(synthX_disc, synthY, '.', 'color', cols(i,:))
end
xlabel('Precision-fixed X', 'FontWeight', 'bold', 'FontSize', 12)
ylabel('Y', 'FontWeight', 'bold', 'FontSize', 12)
set(get(gca, 'Ylabel'), 'Rotation', 0)
xlim([-7, 7])

exportgraphics(gcf,fullfile('figures','simulations_example_gaussian.pdf'),'ContentType','vector')
