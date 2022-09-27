rng('shuffle')
% Moth and muscle to focus on
nmoths = 7;
nmuscles = 10;
knn = 4;
noise = [0, logspace(log10(0.05), log10(6), 120)];
repeats = 150;
n = 5; 

% Preallocate
precision = cell(nmoths, nmuscles);
apriori_precision = cell(nmoths, nmuscles);
MI_disc = cell(nmoths, nmuscles, n);
precision(:) = {zeros(1, n)};
apriori_precision(:) = {zeros(1, n)};
MI_disc(:) = {zeros(length(noise), repeats)};

load('KSG_bin_discretized_data.mat')

% Loop over moths
for m = 6
    % Load data
    load(fullfile('Data',['Moth',num2str(m),'_MIdata.mat']))
    fields = fieldnames(time_data);
    disp(['Moth ',num2str(m)])
    % Loop over muscles
    for muscle = 1:length(fields)
        disp(fields{muscle}(1:end-7))
        % Rename so it's easier to write, get some useful quantities out
        X = time_data.(fields{muscle});
        Y = Tz_WSd;
        
        %---- Bin-discretized real spike data as "fake" X
        % Get bin sizes equivalent to desired max and min noise level
%         minbin = ceil(range(X,'all') / 5);
%         maxbin = round(range(X, 'all') / 1);
%         nbin = round(linspace(minbin, maxbin, n));
        % Preallocate
        fakeX = cell(1,n);
%         bin_edges = cell(1,n);
%         bin = cell(1,n);
        fakeX(:) = {nan(size(X))};
        % Loop over number of bins to discretize by
        for i = 1:n
%             [~, bin_edges{i}, bin{i}] = histcounts(X, nbin(i));
%             fakeX{i} = bin_edges{i}(bin{i} + 1);
%             fakeX{i}(fakeX{i}==bin_edges{i}(1)) = nan;
%             % Indexing above can sometimes cause transpose, so undo if it happens
%             if size(fakeX{i}, 1) < size(fakeX{i}, 2) 
%                 fakeX{i} = fakeX{i}';
%             end
            fakeX{i} = round(X / )
            MI_disc{m,muscle,i} = KSG_precision(fakeX{i}, Y, knn, repeats, noise);
        end
        
        
        % Get precision values
        for i = 1:n
            mis = MI_KSG_subsampling_multispike(fakeX{i}, Y, knn, (1:4));
            mi_sd = findMI_KSG_stddev(mis, size(X,1), false);
            ind = find(mean(MI_disc{m,muscle,i}, 2) < (MI_disc{m,muscle,i}(1,1) - mi_sd), 1);
        %     ind = find(mean(MI_disc{i}, 2) < (MI_disc{i}(1,1) - std(MI_disc{i}(2,:))), 1);
            if isempty(ind)
                precision{m, muscle}(i) = nan;
            else
                precision{m,muscle}(i) = noise(ind);
            end
        end

        apriori_precision{m,muscle} = (max(X, [], 'all') - min(X, [], 'all')) ./ nbin;
    end
end

% % Save results
% save('KSG_bin_discretized_data.mat', 'MI_disc', 'precision', 'apriori_precision', 'noise', 'repeats', 'knn')

% load('KSG_bin_discretized_data.mat')
load('KSG_data.mat')

%%
i = 6;
j = 1;
figure
hold on
cols = copper(n);
for k = 1:n
    meanMI = mean(MI_disc{i,j,k}, 2);
    pad = [nan(1, sg_window), meanMI', nan(1, sg_window)];
    grad = conv(pad, -1 * g(:,2), 'same'); 
    grad = grad(sg_window+1:end-sg_window);
    plot(log10(noise), meanMI, 'color', cols(k,:))
end

%% Get new precision values using other methods
load('KSG_bin_discretized_data.mat')

deriv_thresh_scale = 0.2;
sg_ord = 2; % Savitsky-golay filter order
sg_window = 11; % Savitsky-golay filter window length
[b,g] = sgolay(sg_ord, sg_window);
x = log10(noise);
s = 30; % how many samples on each side to fit line to

deriv_precision = cell(nmoths, nmuscles);
twoline_precision = cell(nmoths, nmuscles);
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
        plot(apriori_precision{i,j}, precision{i,j}, '*', 'color', cols(max_nspike(i,j),:))
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
        plot(apriori_precision{i,j}, deriv_precision{i,j}, '*', 'color', cols(max_nspike(i,j),:))
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
        plot(apriori_precision{i,j}, twoline_precision{i,j}, '*', 'color', cols(max_nspike(i,j),:))
    end
end
plot(get(gca,'xlim'), get(gca,'xlim'), 'k-')
title('twoline')
xlabel('A priori precision (ms)')
ylabel('Measured precision (ms)')

%%

% Fully synthetic dataset 
num_points = 2500;
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


minbin = ceil(range(synthX,'all') / 5);
maxbin = round(range(synthX, 'all') / 1);
nbin = round(linspace(minbin, maxbin, n));

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

%%
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
    mi_sd = findMI_KSG_stddev(mis, num_points, false);
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





