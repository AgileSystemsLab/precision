rng('shuffle')
% Main constants
nmoths = 7;
nmuscles = 10;
knn = 4;
noise = [0, logspace(log10(0.05), log10(6), 120)];
repeats = 150;
prec_levels = 1:0.5:4;
n = length(prec_levels);


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


%% Get precision values, make plots

deriv_thresh_scale = 0.3;
sg_ord = 2; % Savitsky-golay filter order
sg_window = 11; % Savitsky-golay filter window length
[b,g] = sgolay(sg_ord, sg_window);
x = log10(noise);
s = 30; % how many samples on each side to fit line to

precision = cell(nmoths, nmuscles);
deriv_precision = cell(nmoths, nmuscles);
twoline_precision = cell(nmoths, nmuscles);
precision(:) = {nan(1, n)};
deriv_precision(:) = {nan(1, n)};
twoline_precision(:) = {nan(1, n)};
max_nspike = zeros(nmoths, nmuscles);
for i = 1:nmoths
    load(fullfile('Data',['Moth',num2str(i),'_MIdata.mat']))
    fields = fieldnames(time_data);
    for j = 1:nmuscles
        max_nspike(i,j) = size(time_data.(fields{j}), 2);
        for k = 1:n
            meanMI = mean(MI_disc{i,j,k}, 2);
            % Get precision values, OG method
            mis = MI_KSG_subsampling_multispike(fakeX{i,j,k}, Tz_WSd, knn, (1:4));
            mi_sd = findMI_KSG_stddev(mis, size(Tz_WSd,1), false);
            ind = find(mean(MI_disc{i,j,k}, 2) < (MI_disc{i,j,k}(1,1) - mi_sd), 1);
            if isempty(ind)
                precision{i,j}(k) = nan;
            else
                precision{i,j}(k) = noise(ind);
            end
            % Derivative method
            pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
            grad = conv(pad, -1 * g(:,2), 'same'); 
            grad = grad(sg_window+1:end-sg_window);
            deriv_thresh = min(grad) * deriv_thresh_scale;
            deriv_prec_ind = find(grad < deriv_thresh, 1);
            deriv_precision{i,j}(k) = noise(deriv_prec_ind);
            % Two line
            v = pca([x(2:s+1)' meanMI(2:s+1)]);
            beta = v(2,1)/v(1,1);
            int = mean(meanMI(2:s+1)) - beta * mean(x(2:s+1));
            lbeta = [beta, int];
            v = pca([x(end-s:end)' meanMI(end-s:end)]);
            beta = v(2,1)/v(1,1);
            int = mean(meanMI(end-s:end)) - beta * mean(x(end-s:end));
            rbeta = [beta, int];
            twoline_precision{i,j}(k) = 10^((rbeta(2) - lbeta(2)) / (lbeta(1) - rbeta(1)));
        end
    end
end

cols = hot(max(max_nspike, [], 'all'));

figure
hold on
for i = 1:nmoths
    for j = 1:nmuscles
        plot(prec_levels, precision{i,j}, '*', 'color', cols(max_nspike(i,j),:))
    end
end
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
title('original')
xlabel('A priori precision (ms)')
ylabel('Measured precision (ms)')

figure
hold on
for i = 1:nmoths
    for j = 1:nmuscles
        plot(prec_levels, deriv_precision{i,j}, '*', 'color', cols(max_nspike(i,j),:))
    end
end
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
title('derivative')
xlabel('A priori precision (ms)')
ylabel('Measured precision (ms)')

figure
hold on
for i = 1:nmoths
    for j = 1:nmuscles
        plot(prec_levels, twoline_precision{i,j}, '*', 'color', cols(max_nspike(i,j),:))
    end
end
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
title('twoline')
xlabel('A priori precision (ms)')
ylabel('Measured precision (ms)')

%% Fully synthetic dataset from bivariate gaussian

repeats_at_corr = 1;
num_points = 2500;
corr = [0.5, 0.6, 0.7, 0.8, 0.9];
ncorr = length(corr);

% Generate fully synthetic dataset 
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
%%

deriv_thresh_scale = 0.4;

thresh_scale = linspace(0.2, 0.7, 100);
err = zeros(size(thresh_scale));

for ii = 1:length(thresh_scale)
    deriv_thresh_scale = thresh_scale(ii);
    % Plot a priori precision vs actual
    % Find actual precision
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

figure
plot(thresh_scale, err)

% figure
% hold on
% for i = 1:ncorr*repeats_at_corr
%         plot(prec_levels, synth_precision(i,:), '*')
% end
% plot(get(gca, 'xlim'), get(gca,'xlim'), 'k-')
% xlabel('A priori precision')
% ylabel('Observed precision')




