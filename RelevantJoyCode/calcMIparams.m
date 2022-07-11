function MIdata = calcMIparams(spike_data,time_data,torque_data,num_n)

% Code by Joy Putney, last modified 5/22/19

% MIdata = a 1 x 3 matrix containing spike count MI, spike timing MI, and
% the total MI (sum of spike count and spike timing MI.)

% spike_data = a w x 1 matrix containing spike counts in each wing stroke,
% w.

% time_data = a w x S matrix containing spike timings in each wing stroke,
% with a maximum # of spikes, S.

% torque_data = a w x p matrix containing a representation of torque with
% dimensionality p.

% num_n = # of nearest neighbors.

% This sets the dimension of representation of torque_data
num_f = 1:size(torque_data,2);

% Add noise to spike data to prevent weird estimations when spike counts
% are the same for many wing strokes.
spike_data_noised = spike_data + 0.0001*randn(size(spike_data,1),size(spike_data,2));

sMI = MIxnyn(spike_data_noised,torque_data(:,num_f),num_n);
sMI = sMI./log(2); % convert from nats to bits, since C code gives all values in nats.

% Calculate mutual information between timings and flower
spike_counts = unique(spike_data);
for j = 1:length(spike_counts)
    currs = spike_counts(j);
    if currs ~= 0
        idx = find(spike_data == currs);
        prob(j) = length(idx)/size(spike_data,1);
        if length(idx)>num_n+1 && length(idx)>currs && length(idx)>size(torque_data,2)
            KraskovP.tCounts(j) = MIxnyn(time_data(idx,1:currs),torque_data(idx,num_f),num_n);
        else
            KraskovP.tCounts(j) = 0;
        end
    end
end
KraskovP.tCounts = KraskovP.tCounts./log(2); % convert from nats to bits
tCounts = KraskovP.tCounts;
tMI = sum(prob.*KraskovP.tCounts);

totMI = sMI+tMI;

MIdata = [sMI tMI totMI];



