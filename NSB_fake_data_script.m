% NSB fake data script
rng('shuffle')

nmoths = 7;
nmuscles = 10;
nspikebins = 70;
ntorque = 2;
nbiasrepeats = 10;
prec_levels = 1:0.5:4;
n = length(prec_levels);


do_long_runs_real_dataset = false;
do_long_runs_synth_dataset = false;


%% Fixed precision on real dataset
if do_long_runs_real_dataset
    disp('---------- Real dataset ----------')
    % Preallocate
    S_nsbword = cell(nmoths, nmuscles, n);
    dS_nsbword = cell(nmoths, nmuscles, n);
    conditionalentropy = cell(nmoths, nmuscles, n);
    conditionaldS = cell(nmoths, nmuscles, n);
    conditionalentropy_bias = cell(nmoths, nmuscles, n);
    MI_real = cell(nmoths, nmuscles, n);
    MI_real_STD = cell(nmoths, nmuscles, n);
    bias = cell(nmoths, nmuscles, n);
    bins = cell(nmoths, nmuscles, n);
    fakeX = cell(nmoths, nmuscles, n);
    
    S_nsbword(:) = {zeros(1, nspikebins)};
    dS_nsbword(:) = {zeros(1, nspikebins)};
    conditionalentropy(:) = {zeros(1, nspikebins)};
    conditionaldS(:) = {zeros(1, nspikebins)};
    conditionalentropy_bias(:) = {zeros(1, nspikebins)};
    MI_real(:) = {zeros(1, nspikebins)};
    MI_real_STD(:) = {zeros(1, nspikebins)};
    bias(:) = {zeros(1, nspikebins)};
    bins(:) = {zeros(1, nspikebins)};

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
                fakeX{i,j,k} = round(time_data.(fields{j}) / prec_levels(k)) * prec_levels(k);
                % Run main estimator functions
                [S_nsbword{i,j,k}, dS_nsbword{i,j,k}, ~, conditionalentropy{i,j,k}, conditionaldS{i,j,k}, ~] = NSB_precision(fakeX{i,j,k}, Tz_WSd, nspikebins, ntorque);
                [conditionalentropy_bias{i,j,k}, ~, ~, ~] = NSB_precision_bias(fakeX{i,j,k}, Tz_WSd, nspikebins, ntorque, nbiasrepeats);
                % Calculate other quantities
                bins{i,j,k} = range(time_data.(fields{j}), 'all') ./ (1:nspikebins);
                bias{i,j,k} = mean(conditionalentropy_bias{i,j,k}, 1);
                MI_real_STD{i,j,k} = sqrt(dS_nsbword{i,j,k}.^2 + conditionaldS{i,j,k}.^2);
                MI_real{i,j,k} = S_nsbword{i,j,k} - conditionalentropy{i,j,k} - (S_nsbword{i,j,k} - bias{i,j,k});
            end
        end
    end
    % Save results
    save('NSB_sim_real_data_fixed_precision.mat', ...
        'S_nsbword', 'dS_nsbword', 'conditionalentropy', 'conditionaldS', 'conditionalentropy_bias', 'bias', ...
        'MI_real', 'MI_real_STD', 'bins', 'fakeX', 'nspikebins', 'ntorque', 'prec_levels')
end


%% Fixed precision on fully synthetic dataset from bivariate gaussian
repeats_at_corr = 4;
num_points = 2500;
corr = [0.5, 0.6, 0.7, 0.8, 0.9];
ncorr = length(corr);
if do_long_runs_synth_dataset
    disp('---------- Fully synthetic dataset ----------')
    muX = 0;         %X mean of the bivariate gaussian
    muY = 0;         %Y mean of the bivariate gaussian
    sigmaX = 2;      %X standard deviation   
    sigmaY = 2;      %Y standard deviation
    
    % preallocate
    S_nsbword = cell(ncorr, repeats_at_corr, n);
    dS_nsbword = cell(ncorr, repeats_at_corr, n);
    conditionalentropy = cell(ncorr, repeats_at_corr, n);
    conditionaldS = cell(ncorr, repeats_at_corr, n);
    conditionalentropy_bias = cell(ncorr, repeats_at_corr, n);
    MI_synth = cell(ncorr, repeats_at_corr, n);
    MI_synth_STD = cell(ncorr, repeats_at_corr, n);
    bias = cell(ncorr, repeats_at_corr, n);
    bins = cell(ncorr, repeats_at_corr, n);
    
    S_nsbword(:) = {zeros(1, nspikebins)};
    dS_nsbword(:) = {zeros(1, nspikebins)};
    conditionalentropy(:) = {zeros(1, nspikebins)};
    conditionaldS(:) = {zeros(1, nspikebins)};
    conditionalentropy_bias(:) = {zeros(1, nspikebins)};
    MI_real(:) = {zeros(1, nspikebins)};
    MI_real_STD(:) = {zeros(1, nspikebins)};
    bias(:) = {zeros(1, nspikebins)};
    bins(:) = {zeros(1, nspikebins)};
    for i = 1:ncorr
        for j = 1:repeats_at_corr
            %generate a correlated X and Y
            x1 = normrnd(0, 1, num_points, 1);
            x2 = normrnd(0, 1, num_points, 1);
            x3 = corr(i) .* x1 + (1 - corr(i)^2)^.5 .* x2;
            synthX = muX + x1 * sigmaX;
            synthY = muY + x3 * sigmaY;
            % Loop over fixed precision levels
            for k = 1:n
                synthX_disc = round(synthX / prec_levels(k)) * prec_levels(k);
                % Run main estimator functions
                [S_nsbword{i,j,k}, dS_nsbword{i,j,k}, ~, conditionalentropy{i,j,k}, conditionaldS{i,j,k}, ~] = NSB_precision(synthX_disc, synthY, nspikebins, ntorque);
                [conditionalentropy_bias{i,j,k}, ~, ~, ~] = NSB_precision_bias(synthX_disc, synthY, nspikebins, ntorque, nbiasrepeats);
                % Calculate other quantities
                bins{i,j,k} = range(synthX_disc, 'all') ./ (1:nspikebins);
                bias{i,j,k} = mean(conditionalentropy_bias{i,j,k}, 1);
                MI_synth_STD{i,j,k} = sqrt(dS_nsbword{i,j,k}.^2 + conditionaldS{i,j,k}.^2);
                MI_synth{i,j,k} = S_nsbword{i,j,k} - conditionalentropy{i,j,k} - (S_nsbword{i,j,k} - bias{i,j,k});
            end
        end
    end
    
    % Save results
    save('NSB_sim_synthetic_data_fixed_precision.mat', ...
        'S_nsbword', 'dSnsbword', 'conditionalentropy', 'conditionaldS', 'conditionalentropy_bias', ...
        'MI_synth', 'MI_synth_STD', 'bias', 'bins', ...
        'corr', 'repeats_at_corr', 'num_points')
end


%% Methods comparison
load('NSB_sim_real_data_fixed_precision.mat')
% load('NSB_sim_synthetic_data_fixed_precision.mat')


sg_ord = 2; % Savitsky-golay filter order
sg_window = 13; % Savitsky-golay filter window length
[b,g] = sgolay(sg_ord, sg_window);
min_deriv_thresh = 0.2;
deriv_thresh = 0.04;

cols = copper(n);
pi = 5; pj = 9;

precision_max = nan(nmoths, nmuscles, n);
precision_deriv = nan(nmoths, nmuscles, n);
precision_other = nan(nmoths, nmuscles, n);
for i = 1:nmoths
    for j = 1:nmuscles
        if (i == pi) && (j == pj)
            figure
            ax1 = subplot(2, 1, 1);
            hold on
            ax2 = subplot(2, 1, 2);
            hold on
        end
        for k = 1:n
            % Max value method
            useMI = S_nsbword{i,j,k} - conditionalentropy{i,j,k};
            ind = find(useMI(1:end-1) == useMI(2:end), 1);
            precision_max(i,j,k) = bins{i,j,k}(ind);
%             [~,ind] = max(S_nsbword{i,j,k} - conditionalentropy{i,j,k});
%             precision_max(i,j,k) = bins{i,j,k}(ind);
            % Derivative zero-crossing method
            pad = [nan(1, sg_window), MI_real{i,j,k}, nan(1, sg_window)];
            grad = conv(pad, 1 * g(:,2), 'same'); 
            grad = grad(sg_window+1:end-sg_window);
            inds = find((sign(grad(2:end))==-1) & (sign(grad(1:end-1))==1));
%             inds = find(grad/min(grad) > min_deriv_thresh);
            if isempty(inds)
                precision_deriv(i,j,k) = nan;
            else
                precision_deriv(i,j,k) = bins{i,j,k}(inds(1));
            end
            % Other method
            useMI = MI_real{i,j,k};
            mwind = zeros(1, length(MI_real{i,j,k}));
            for ii = 1:length(MI_real{i,j,k})
%                 mwind(ii) = mean(MI_real{i,j,k}(ii:end));
                mwind(ii) = mean(useMI(ii:end));
            end
%             [~,ind] = max(MI_real{i,j,k}-mwind);
            [~,ind] = max(useMI-mwind);
            precision_other(i,j,k) = bins{i,j,k}(ind);

            if (i == pi) && (j == pj)
                plot(ax1, bins{i,j,k}, S_nsbword{i,j,k} - conditionalentropy{i,j,k}, 'color', cols(k,:))
%                 plot(ax1, bins{i,j,k}(ind), useMI(ind), 'r.')
                plot(ax2, bins{i,j,k}, MI_real{i,j,k}, 'color', cols(k,:))
            end
        end
    end
end
set(ax1, 'Xscale', 'log')
set(ax2, 'Xscale', 'log')
% %%
% i=2;j=2;k=5;
% mwind = zeros(1, length(MI_real{i,j,k}));
% for ii = 1:length(MI_real{i,j,k})
%     mwind(ii) = mean(MI_real{i,j,k}(ii:end));
% end
% 
% figure
% hold on
% % plot(MI_real{i,j,k})
% plot(MI_real{i,j,k}-mwind)
% % plot(mwind)

% Reshape for plotting
precision_max = reshape(precision_max, [nmoths*nmuscles, n]);
precision_deriv = reshape(precision_deriv, [nmoths*nmuscles, n]);
precision_other = reshape(precision_other, [nmoths*nmuscles, n]);
% Plot
figure('outerposition', [440 366 967 300])
t = tiledlayout(1, 3);
col = '#4472C4';
% Max value method
nexttile()
hold on
errorbar(prec_levels, mean(precision_max, 1, 'omitnan'), std(precision_max, 1, 'omitnan'), 'o', ...
    'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
    'LineWidth', 1, 'CapSize', 13)
xlim([prec_levels(1)-0.25, prec_levels(end)+0.25])
ylim([0, 5])
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
title('Max Value Method')
xlabel('Actual precision (ms)')
ylabel('Measured precision (ms)')
% Derivative zero-crossing method
nexttile()
hold on
errorbar(prec_levels, mean(precision_deriv, 1, 'omitnan'), std(precision_deriv, 1, 'omitnan'), 'o', ...
    'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
    'LineWidth', 1, 'CapSize', 13)
xlim([prec_levels(1)-0.25, prec_levels(end)+0.25])
ylim([0, 5])
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
title('Deriv zero-cross Method')
xlabel('Actual precision (ms)')
ylabel('Measured precision (ms)')
% other method
nexttile()
hold on
errorbar(prec_levels, mean(precision_other, 1, 'omitnan'), std(precision_other, 1, 'omitnan'), 'o', ...
    'color', col, 'Marker', 'o', 'MarkerEdgeColor', col, 'MarkerFaceColor', col, ...
    'LineWidth', 1, 'CapSize', 13)
xlim([prec_levels(1)-0.25, prec_levels(end)+0.25])
ylim([0, 5])
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
title('Other Method')
xlabel('Actual precision (ms)')
ylabel('Measured precision (ms)')


%% Distributions of precision at specific known levels
% Actually find precision for both real and synthetic fixed datasets
load('NSB_sim_real_data_fixed_precision.mat')
load('NSB_sim_synthetic_data_fixed_precision.mat')
min_peak_height = 0.2;
min_peak_dist = 5;

sg_ord = 2; % Savitsky-golay filter order
sg_window = 13; % Savitsky-golay filter window length
[b,g] = sgolay(sg_ord, sg_window);
% Real
figure
hold on
cols = copper(n);
precision_real = zeros(nmoths, nmuscles, n);
for i = 1:nmoths
    for j = 1:nmuscles
        for k = 1:n
            pad = [nan(1, sg_window), MI_real{i,j,k}, nan(1, sg_window)];
            grad = conv(pad, -1 * g(:,2), 'same'); 
            grad = grad(sg_window+1:end-sg_window);
            dpad = [nan(1, sg_window), grad, nan(1, sg_window)];
            dgrad = conv(dpad, -1 * g(:,2), 'same'); 
            dgrad = dgrad(sg_window+1:end-sg_window);
            [~,inds] = findpeaks(dgrad/min(dgrad), 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_dist);
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
            pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
            grad = conv(pad, -1 * g(:,2), 'same'); 
            grad = grad(sg_window+1:end-sg_window);
            dpad = [nan(1, sg_window), grad, nan(1, sg_window)];
            dgrad = conv(dpad, -1 * g(:,2), 'same'); 
            dgrad = dgrad(sg_window+1:end-sg_window);
            [~,inds] = findpeaks(dgrad/min(dgrad), 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_dist);
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
