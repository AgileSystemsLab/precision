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

synth_repeats = 50;
num_points = 2500;

% Generate fully synthetic dataset 
muX = 0;         %X mean of the bivariate gaussian
muY = 0;         %Y mean of the bivariate gaussian
sigmaX = 2;      %X standard deviation   
sigmaY = 2;      %Y standard deviation
correlation = 0.9;

% preallocate
synthX_disc = cell(synth_repeats,n);
synthY = cell(synth_repeats, 1);
MI_synth = cell(synth_repeats,n);
synthX_disc(:) = {nan(num_points, 1)};
synthY(:) = {nan(num_points, 1)};
MI_synth(:) = {zeros(length(noise), repeats)};
for r = 1:synth_repeats
    %generate a correlated X and Y
    x1 = normrnd(0, 1, num_points, 1);
    x2 = normrnd(0, 1, num_points, 1);
    x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
    synthX = muX + x1 * sigmaX;
    synthY{r} = muY + x3 * sigmaY;
    % Loop over fixed precision levels
    for i = 1:n
        synthX_disc{r,i} = round(synthX / prec_levels(i)) * prec_levels(i);
        MI_synth{r,i} = KSG_precision(synthX_disc{r,i}, synthY{r}, knn, repeats, noise);
    end
end

% Save results
save('KSG_sim_synthetic_data_fixed_precision.mat', 'MI_synth', 'synthX_disc', 'synthY', ...
    'noise', 'repeats', 'knn', 'synth_repeats', 'num_points', 'correlation')

%%

i = 3;

figure
hold on
cols = copper(n);
for j = 1:n
    plot(log10(noise), mean(MI_synth{i,j}, 2), 'color', cols(j,:))
end

%%
% % Plot MI vs noise curves
% figure
% hold on
% cols = copper(n);
% for i = 1:n
%     plot(log10(noise), mean(MI_synth{i}, 2), 'color', cols(i,:));
% end
% MI_synth_real = KSG_precision(synthX, synthY, knn, repeats, noise);
% mseb(log10(noise), mean(MI_synth_real,2), std(MI_synth_real, 0, 2)');
% xlabel('log10(noise)')
% ylabel('MI (bits)')

% Plot a priori precision vs actual
synth_precision = zeros(synth_repeats, n);
synth_deriv_precision = zeros(synth_repeats, n);
for r = 1:synth_repeats
    for i = 1:n
        meanMI = mean(MI_synth{r,i}, 2);
        % Original method
        mis = MI_KSG_subsampling_multispike(synthX_disc{r,i}, synthY{r}, knn, (1:4));
        mi_sd = findMI_KSG_stddev(mis, num_points, false);
        ind = find(meanMI < (MI_synth{r,i}(1,1) - mi_sd), 1);
        synth_precision(r,i) = noise(ind);
        % Derivative method
        pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
        grad = conv(pad, -1 * g(:,2), 'same'); 
        grad = grad(sg_window+1:end-sg_window);
        deriv_thresh = min(grad) * deriv_thresh_scale;
        ind = find(grad < deriv_thresh, 1);
        synth_deriv_precision(r,i) = noise(ind);
    end
end

figure
hold on
for i = 1:synth_repeats
    plot(prec_levels, synth_precision(i,:), '*')
end
plot(get(gca, 'xlim'), get(gca,'xlim'), 'k-')
title('original')
xlabel('A priori precision')
ylabel('Observed precision')
figure
hold on
for i = 1:synth_repeats
    plot(prec_levels, synth_deriv_precision(i,:), '*')
end
plot(get(gca, 'xlim'), get(gca,'xlim'), 'k-')
title('derivative')
xlabel('A priori precision')
ylabel('Observed precision')


% %% Fully synthetic dataset, vary correlation factor
% corr = [0.4, 0.6, 0.8, 0.9];
% 
% num_points = 2700;
% muX = 0;         %X mean of the bivariate gaussian
% muY = 0;         %Y mean of the bivariate gaussian
% sigmaX = 2;      %X standard deviation   
% sigmaY = 2;      %Y standard deviation
% %generate a correlated X and Y
% x1 = normrnd(0, 1, num_points, 1);
% x2 = normrnd(0, 1, num_points, 1);
% 
% figure(1)
% hold on
% figure(2)
% hold on
% cols = copper(length(corr));
% 
% 
% 
% for j = 1:length(corr)
%     correlation = corr(i);
%     x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
%     synthX = muX + x1 * sigmaX;
%     synthY = muY + x3 * sigmaY;
%     MI = KSG_precision(synthX, synthY, knn, repeats, noise);
%     figure(1)
%     mseb(log10(noise), mean(MI, 2), std(MI, 0, 2)', struct('col', {{cols(i,:)}}));
%     synthX_disc = cell(1,n);
%     MI_synth = cell(1,n);
%     bin_edges = cell(1,n);
%     bin = cell(1,n);
%     synthX_disc(:) = {nan(size(synthX))};
%     MI_synth(:) = {zeros(length(noise), repeats)};
%     % Loop over bins to discretize by
%     for i = 1:n
%         [~, bin_edges{i}, bin{i}] = histcounts(synthX, nbin(i));
%         synthX_disc{i} = bin_edges{i}(bin{i} + 1);
%         synthX_disc{i}(synthX_disc{i}==bin_edges{i}(1)) = nan;
%         % Indexing above can sometimes cause transpose, so undo if it happens
%         if size(synthX_disc{i}, 1) < size(synthX_disc{i}, 2) 
%             synthX_disc{i} = synthX_disc{i}';
%         end
%         MI_synth{i} = KSG_precision(synthX_disc{i}, synthY, knn, repeats, noise);
%     end
% 
%     figure(2)
% 
% 
% end





