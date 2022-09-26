rng('shuffle')
% Moth and muscle to focus on
moth = '3';
muscle = 'LDLM';
knn = 4;
noise = [0, logspace(log10(0.05), log10(6), 120)];
repeats = 150;
n = 10; 

%---- Load data
load(fullfile('Data',['Moth',num2str(moth),'_MIdata.mat']))
% Rename so it's easier to write, get some useful quantities out
X = time_data.([muscle,'strokes']);
Y = Tz_WSd;
Nspike = sum(~isnan(X), 2);
unq = unique(Nspike);

%% Real data bin-discretized 
%---- Bin-discretized real spike data as "fake" X
% Get bin sizes equivalent to desired max and min noise level
minbin = ceil(range(X,'all') / 5);
maxbin = round(range(X, 'all') / 1);
nbin = round(linspace(minbin, maxbin, n));
% Preallocate
MI_disc = cell(1,n);
fakeX = cell(1,n);
bin_edges = cell(1,n);
bin = cell(1,n);
MI_disc(:) = {zeros(length(noise), repeats)};
fakeX(:) = {nan(size(X))};
% Loop over number of bins to discretize by
for i = 1:n
    [~, bin_edges{i}, bin{i}] = histcounts(X, nbin(i));
    fakeX{i} = bin_edges{i}(bin{i} + 1);
    fakeX{i}(fakeX{i}==bin_edges{i}(1)) = nan;
    % Indexing above can sometimes cause transpose, so undo if it happens
    if size(fakeX{i}, 1) < size(fakeX{i}, 2) 
        fakeX{i} = fakeX{i}';
    end
    MI_disc{i} = KSG_precision(fakeX{i}, Y, knn, repeats, noise);
end
MI_real = KSG_precision(X, Y, knn, repeats, noise);


% Get precision values
precision = zeros(1, n);
precision_ind = zeros(1, n);
for i = 1:n
    mis = MI_KSG_subsampling_multispike(fakeX{i}, Y, knn, (1:4));
    mi_sd = findMI_KSG_stddev(mis, size(X,1), false);
    precision_ind(i) = find(mean(MI_disc{i}, 2) < (MI_disc{i}(1,1) - mi_sd), 1);
%     precision_ind(i) = find(mean(MI_disc{i}, 2) < (MI_disc{i}(1,1) - std(MI_disc{i}(2,:))), 1);
    precision(i) = noise(precision_ind(i));
end

% Plot a priori precision vs observed
apriori_precision = (max(X, [], 'all') - min(X, [], 'all')) ./ nbin;
figure
hold on
plot(apriori_precision, precision, '*')
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')

% Plot MI against noise at different binning levels
figure
hold on
cols = copper(n);
for i = 1:n
    % Plot line of MI against noise
    plot(log10(noise), mean(MI_disc{i}, 2), 'color', cols(i,:))
    % Plot precision points
%     plot(log10(noise(precision_ind(i))), mean(MI_disc{i}(precision_ind(i),:)), '*')
%     plot(log10(noise(deriv_prec_ind(i))), mean(MI_disc{i}(deriv_prec_ind(i),:)), 'o')
end
mseb(log10(noise), mean(MI_real,2), std(MI_real, 0, 2)');
xlabel('log10(noise)')
ylabel('MI (bits)')
title(['Moth ', moth, ' ', muscle])


%%

% Fully synthetic dataset 
num_points = 2700;
muX = 0;         %X mean of the bivariate gaussian
muY = 0;         %Y mean of the bivariate gaussian
sigmaX = 2;      %X standard deviation   
sigmaY = 2;      %Y standard deviation
correlation = 0.9;
%generate a correlated X and Y
x1 = normrnd(0, 1, num_points, 1);
x2 = normrnd(0, 1, num_points, 1);
x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
synthX = muX + x1 * sigmaX;
synthY = muY + x3 * sigmaY;

synthX_disc = cell(1,n);
MI_synth = cell(1,n);
bin_edges = cell(1,n);
bin = cell(1,n);
synthX_disc(:) = {nan(size(synthX))};
MI_synth(:) = {zeros(length(noise), repeats)};
% Loop over bins to discretize by
for i = 1:n
    [~, bin_edges{i}, bin{i}] = histcounts(synthX, nbin(i));
    synthX_disc{i} = bin_edges{i}(bin{i} + 1);
    synthX_disc{i}(synthX_disc{i}==bin_edges{i}(1)) = nan;
    % Indexing above can sometimes cause transpose, so undo if it happens
    if size(synthX_disc{i}, 1) < size(synthX_disc{i}, 2) 
        synthX_disc{i} = synthX_disc{i}';
    end
    MI_synth{i} = KSG_precision(synthX_disc{i}, synthY, knn, repeats, noise);
end
MI_synth_real = KSG_precision(synthX, synthY, knn, repeats, noise);


% Plot MI vs noise curves
figure
hold on
cols = copper(n);
for i = 1:n
    plot(log10(noise), mean(MI_synth{i}, 2), 'color', cols(i,:));
end
mseb(log10(noise), mean(MI_synth_real,2), std(MI_synth_real, 0, 2)');
xlabel('log10(noise)')
ylabel('MI (bits)')

% Plot a priori precision vs actual
synth_precision = zeros(1, n);
apriori_precision = (max(synthX, [], 'all') - min(synthX, [], 'all')) ./ nbin;
for i = 1:n
    mis = MI_KSG_subsampling_multispike(synthX_disc{i}, synthY, knn, (1:4));
    mi_sd = findMI_KSG_stddev(mis, size(X,1), false);
    ind = find(mean(MI_synth{i}, 2) < (MI_synth{i}(1,1) - mi_sd), 1);
    synth_precision(i) = noise(ind);
end

figure
hold on
plot(apriori_precision, synth_precision, '*')
plot(get(gca, 'xlim'), get(gca,'xlim'), 'k-')
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





