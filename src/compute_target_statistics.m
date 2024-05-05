function [true_statistics] = compute_target_statistics(target_train_name,obj_configs)
%% Compute the activity statistics of the target spike trains.
%   -Input
%     obj_configs: string; filename of the spike trains to customize network model to.  
%     obj_configs: struct, configurations for customization.

%   -Output
%      true_statistics: table, the six activity statistics of the target data


load(strcat('./data/',target_train_name,'.mat'));
n_sampling = obj_configs.n_sampling;
Tw=obj_configs.Tw;
Tburn=obj_configs.Tburn;
Ne1=obj_configs.Ne1;
n_neuron = 50; %fixed for the short demo
dim_method = obj_configs.dim_method;
Ic1 = sample_e_neurons(spike_train,Ne1,1);
[rate1,var1, FanoFactor1, mean_corr1, unstable_flag, sampling_inds, re,low_rate_flag]=compute_stats(spike_train,Ic1,Tw,Tburn,n_sampling, n_neuron, 1);
if low_rate_flag
    error('low firing rate, use a different spike train')
end
if unstable_flag
    error('unrealistic spiking, use a different spike train')
end


[fa_percentshared, fa_normevals, fa_dshared] = compute_pop_stats(sampling_inds, re, n_neuron, Tw, dim_method);


% In the manuscript, we compute the variance of each statistic over
% multiple simulations of the same network parameter set. Here we preset
% the variance using empirical numbers for simplicity.

rate_mean=rate1;
rate_var=0.1;

fano_mean=FanoFactor1; 
fano_var=0.0002;             

mean_corr_mean=mean_corr1; 
mean_corr_var=2e-6;    

fa_percent_mean=fa_percentshared;
fa_percent_var=1e-6;    

fa_dim_mean=fa_dshared; 
fa_dim_var=0.1;   

fa_normeval_mean=fa_normevals;  
fa_normeval_var=0.01;

default_weights=ones(1,6);
true_statistics=table(n_neuron,rate_mean,rate_var,fano_mean,fano_var,mean_corr_mean,mean_corr_var,fa_percent_mean,fa_percent_var,fa_dim_mean,fa_dim_var,fa_normeval_mean,fa_normeval_var,default_weights);
true_statistics
%output_file_name = strcat(target_train_name,'_stats');
%save(strcat('./data/',output_file_name,'.mat'),'true_statistics');
end 