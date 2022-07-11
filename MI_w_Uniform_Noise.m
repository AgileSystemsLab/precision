function [miVec, miextra] = MI_w_Uniform_Noise(motoroutput, spikedata, spikerate)

%% General Information

% To run this code, you will also need the KSG MI Estimator package
% available at this repository: 
% https://github.com/EmoryUniversityTheoreticalBiophysics/ContinuousMIEstimation

% When using this code, please cite the following:

% Putney, J., Niebur, N., Conn, R., & Sponberg, S. An information theoretic
% method to resolve millisecond-scale spike timing precision in a 
% comprehensive motor program. bioRxiv 2021.07.14.452403.

% Holmes, CM. & Nemenman, I. (2019). Estimation of mutual information for
% real-valued data with error bars and controlled bias. Phys Rev E, 100(2), 
% 022404.

% Srivastava, K., Holmes, CM., Vellema, M., Pack, AR., Coen, PHE.,
% Nemenman, I., & Sober, SJ. (2017). Motor control by precisely timed spike
% patterns. PNAS 114(5), 1171-1176.

% Kraskov, A., Stoegbauer, H., & Grassberger, P. (2004). Estimating
% mutual information. Phys Rev E, 69(6), 066138.

%% Variables

% motoroutput - the motor output, with N x L dimensions where N is the number
% of observations and L is the number of variables used to describe the
% motor output. In Putney et al., L is 2 for a 2-dim PCA scores
% representation. This variable can be anything the spiking activity is
% encoding, including information about a sensory stimulus.

% spikedata- the set of spikes with timing data included, with N x C_max dimesnions.
% This should be a matrix, where the rows correspond to each observation,
% and the columns denoting the time of the spikes in ms. If variable #s of
% spikes occur in each observation, NaNs can be used for observations where
% fewer than C_max spikes occurs. In Putney et al., each observation is a
% wing stroke.

% spikerate - This should be a N x 1 matrix, where the rows correspond to
% each observation and the entires denote the number of spikes in that
% observation.

% miVec- A 150 x length(noise) matrix of the mutual information (MI) in
% bits/observation. Here, the MI is estimated at each value of uniform noise
% in the vector noisevec 150 times. The output of the function is therefore
% the values of MI at each level of noise defined in noisevec, where the
% rows are 150 runs of that value.

numofneighbors = 4; % number of k nearest neighbors for KSG MI estimator
noise = (0:.05:6); % range of noise distribution added to the timing data in ms
miVec = zeros(150, length(noise));% preallocating the mivec output, 150 trials at each level of noise

%% Estimating MI with Uniform Noise

minnumspike = min(spikerate); %smallest number of spikes in this dataset
maxnumspike = max(spikerate); %largest number of spikes in this dataset
[j,k] = size(spikedata);

% iterate from smallest to largest number of spikes to create a distribution
% of spikes by the number of spikes, basically separating the torque and
% timing by the number of spikes
for i = minnumspike:maxnumspike
    if i == 0
        spikestructnorm.(['spike', num2str(i)]).timedata = spikedata(spikerate == i, :);
    else
        if k<=i
            spikestructnorm.(['spike', num2str(i)]).timedata = spikedata(spikerate ==  i, 1:k);
        else
            spikestructnorm.(['spike', num2str(i)]).timedata = spikedata(spikerate ==  i, 1:i);
        end
    end
    
    spikestructnorm.(['spike', num2str(i)]).encodedvar = motoroutput(spikerate == i, :);
    spikestructnorm.(['spike', num2str(i)]).probability = length(spikestructnorm.(['spike', num2str(i)]).timedata)/length(spikerate);
    
end

[a,b] = size(miVec);
fn =  fieldnames(spikestructnorm);
miextra = zeros(150, length(noise), length(fn));
% Apply KSG MI estimator to the structure created above
for x = 1:b
    for y = 1:a
        for ffn = fn'
            ffn;
            dbstop if error
            numwordsnoised = spikestructnorm.(ffn{1}).timedata + (noise(x).*rand(size(spikestructnorm.(ffn{1}).timedata,1),size(spikestructnorm.(ffn{1}).timedata,2)));
            [j,k] = size(numwordsnoised);
            if j > numofneighbors && j+1 > str2num(ffn{1}(end))
                % This line calls MIxnyn which is included in the KSG
                % GITHUB repository cited above from Holmes et al. 2019
                wMI = MIxnyn_matlab(numwordsnoised, spikestructnorm.(ffn{1}).encodedvar, numofneighbors, pwd);
                wMI = wMI./log(2); % convert from nats to bits, since C code gives all values in nats.
            else
                
                wMI = 0;
                
            end
            miVec(y,x) = wMI*spikestructnorm.(ffn{1}).probability+miVec(y,x);
            miextra(y,x,str2double(ffn{1}(end))+1) = wMI * spikestructnorm.(ffn{1}).probability;
        end
    end
end

end