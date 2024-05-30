function [rate0, var0, FanoFactor0, mean_corr0, unstable_flag, sampling_inds, re, low_rate_flag] = compute_stats(s1, Ic1, Tw, Tburn, n_sampling, n_neuron, check_stability) 
%% Wrapper function for computing the single-neuron and pairwise statistics of the given spike trains
%   -Input
%     s1: [number spikes, 2], spike trains
%     Ic1: [number of E nerons], indices of the E neurons in the recurrent layer
%     Tw: int, spike count window size
%     Tburn: int, length of spike trains to be burnt (for remove transcients)
%     n_sampling: int, number of neuron samplings for averaging the activity stats
%     n_neuron: int, number of neurons per sampling
%     check_stability: {0, 1}, check if the spiking activity is stable (see Feasibility constraints in Methods)

%   -Output
%      rate0: float;  mean firing rate
%      var0: float, mean variance
%      FanoFactor0: float, mean Fanofactor
%      mean_corr0: float, mean rsc
%      unstable_flag: {0, 1}, if the spiking activity is stable
%      sampling_inds: [number of samplings, n_neuron], neuron sampling indices
%      re: [number of neurons, number of bins], spike count matrix
%      low_rate_flag: {0, 1}, if the spike trains have extreme firing rate



low_rate_flag=0;
unstable_flag=0;

rate_th = 0.5; % low fr threshold
s1=s1(:, ismember(s1(2,:), Ic1));
T_total=max(s1(1,:));
s1=s1(:,s1(1,:)>=Tburn&s1(1,:)<=T_total); 
s1(1,:)=s1(1,:)-Tburn;

re=spktime2count(s1, unique(s1(2,:)), Tw, floor((T_total-Tburn)/Tw), 0); % bin the spike trains into the spike count matrix
% Check feasibility constraint: long transcient
unstable_flag=is_unstable(re);
unstable_flag=0;%%!!!
if unstable_flag
    rate0=NaN;
    var0=NaN;
    FanoFactor0=NaN;
    mean_corr0=NaN;
    sampling_inds=NaN;
    re=NaN;
    warning('long transient oscillation, unstable simulation!')
    return 
end

re=re(mean(re,2)*1000/Tw>rate_th, :);
if check_stability
    % Check feasibility constraint: low fr
    low_rate_flag= (size(re,1)<n_neuron) | (mean(mean(re,2))*1000/Tw>60) ; %mar 27: add high thres rejection
    if low_rate_flag
        rate0=NaN;
        var0=NaN;
        FanoFactor0=NaN;
        mean_corr0=NaN;
        sampling_inds=NaN;
        re=NaN;
        warning('extreme firing rate!')
        return 
    end
end

% Randomly sample n_neuron for n_sampling times to compute the activity statistic to mitigate finite network size var
[sampling_inds]=random_sample_neurons(1:size(re,1), n_sampling, n_neuron);
rate0s=zeros(n_sampling, 1);
var0s=zeros(n_sampling, 1);
FanoFactor0s=zeros(n_sampling, 1);
mean_corr0s=zeros(n_sampling, 1);
for cnt = 1:n_sampling
    tmp = re(sampling_inds(cnt,:), :);
    [rate0,var0, FanoFactor0, mean_corr0]=compute_statistics_only(tmp);
    rate0s(cnt)=rate0*1000/Tw; % look at all neurons for rate
    var0s(cnt)=var0;
    FanoFactor0s(cnt)=FanoFactor0;
    mean_corr0s(cnt)=mean_corr0;
end

rate0=mean(rate0s, 1);
var0=mean(var0s, 1);
FanoFactor0=mean(FanoFactor0s, 1);
mean_corr0=mean(mean_corr0s, 1);

if isnan(rate0) || isnan(var0) || isnan(FanoFactor0) || isnan(mean_corr0)
    error('Invalid stats, possibly due to low firing rate and <=1 neuron passes the rate threshold!')
end

end
