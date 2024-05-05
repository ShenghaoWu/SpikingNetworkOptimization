function [sampling_inds]=random_sample_neurons(Ic1,n_sampling, n_neuron) 
%% Function for randomly sampling neurons from the indices (Ic1), returns the sampling indices.
%   -Input
%      Ics: [number of neurons], neuron indices for sampling
%      n_sampling: int, number of samplings.
%      n_neuron: int, number of neurons per sampling

%   -Output
%      sampling_inds: [n_sampling, n_neuron], the sampling indices.


sampling_inds = zeros(n_sampling, n_neuron);
for cnt = 1:n_sampling
    sampling_inds(cnt,:) = randsample(Ic1,n_neuron);
end

end