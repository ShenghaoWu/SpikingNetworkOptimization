function [objective] = cost_func_short_demo(x, is_surrogate, obj_configs, simulator)
%% Cost function for demo_short.m. Note that this function is similar to cost_func.m, with simplified functionality for demonstration purpose. 
%   -Input
%      x: table, parameter sets for the simulation
%      is_surrogate: {0,1}, if use a short simulation to determine
%      feasibility
%      obj_configs: struct, configurations for the network simulation
%      simulator: function, wrapper for the network simulation, corresponds to spatial_nn_simulation_weight.m in cost_func.m        

%   -Output
%      objective: float; cost of x


%initialization

rng("default");
tic

T=obj_configs.T;
T0=ceil(obj_configs.T/10);
n_sampling = obj_configs.n_sampling;
Tw=obj_configs.Tw;
Tburn=obj_configs.Tburn;
Ne1=obj_configs.Ne1;
true_statistics=obj_configs.true_statistics;
tolerance=obj_configs.tolerance;
stats_filename=obj_configs.stats_filename;
save_stats=obj_configs.save_stats;
statistics_group=obj_configs.statistics_group;
try 
	n_neuron = obj_configs.n_neuron;
catch
	n_neuron=true_statistics.n_neuron;
end
dim_method = obj_configs.dim_method;

%%%%%%%%%%%%%%%%%%%%%%%%
%%% surrogate stats %%%
%%%%%%%%%%%%%%%%%%%%%%%%


warning_flag=0;
if is_surrogate
	try
		[s1] =  simulator(x,T0);
		% if the simulation is terminated not too earlier from T, we still compute its stats
		if max(s1(1,:))<0.9*T0
			error('simulation terminated due to max spikes, parameters discarded')
		end
		Ic1 = sample_e_neurons(s1,Ne1,1);
		[rate0,var0, FanoFactor0, ~, ~, ~, ~,low_rate_flag]=compute_stats(s1,Ic1,Tw,Tburn,n_sampling,n_neuron,0);
		if low_rate_flag
			error('low firing rate in surrogate simulation!')
        end
		 
        cost_surrogate=mean([true_statistics.default_weights(1)*(rate0-true_statistics.rate_mean)^2/true_statistics.rate_var,...
                            true_statistics.default_weights(2)*(FanoFactor0-true_statistics.fano_mean)^2/true_statistics.fano_var]);

		if isnan(cost_surrogate)||cost_surrogate>tolerance
			warning_flag=1;
		end
	catch 
		warning_flag=1;
		warning('something wrong happended during surrogate stage')
	end
else
	cost_surrogate=NaN;
	rate0=NaN;
	var0=NaN;
	FanoFactor0=NaN;
end
%%%%%saving stats (even NaNs)
objective=NaN;
if save_stats
	surrogate_time=toc;
	full_stats_time=0;
	execution_time=table(surrogate_time,full_stats_time);
	[rate1,var1,FanoFactor1,mean_corr1,fa_percentshared,fa_dshared,fa_normevals] = nan_stats(n_neuron);
	full_stats_to_record = generate_full_stats_to_record(statistics_group,objective,rate1,var1,FanoFactor1,mean_corr1,fa_percentshared,fa_dshared,fa_normevals);
	if warning_flag
		try
			surrogate_stats=table(cost_surrogate,rate0,var0,FanoFactor0);
		catch
			cost_surrogate=NaN;
			rate0=NaN;
			var0=NaN;
			FanoFactor0=NaN;
			surrogate_stats=table(cost_surrogate,rate0,var0,FanoFactor0);
		end

		stats_weights=true_statistics.default_weights;
		save_statistics(stats_weights,x,full_stats_to_record,surrogate_stats,execution_time,stats_filename);
		return
	end
end 

if warning_flag
	return
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%% full stats %%%
%%%%%%%%%%%%%%%%%%%%%%%%
tic

warning_flag=0; 

try
	[s1] =  simulator(x,T);
	if max(s1(1,:))<0.9*T
		error('simulation terminated due to max spikes, parameters discarded')
	end
	Ic1 = sample_e_neurons(s1,Ne1,1);
	[rate1,var1, FanoFactor1, mean_corr1, unstable_flag, sampling_inds, re,low_rate_flag]=compute_stats(s1,Ic1,Tw,Tburn,n_sampling, n_neuron, 1);
	if low_rate_flag
		error('low firing rate in main simulation!')
	end

	if unstable_flag
		
		[s1] =  simulator(x,T);
		if max(s1(1,:))<0.9*T
			error('simulation terminated due to max spikes, parameters discarded')
		end
		Ic1 = sample_e_neurons(s1,Ne1,1);
		[rate1,var1, FanoFactor1, mean_corr1, unstable_flag, sampling_inds, re, low_rate_flag]=compute_stats(s1,Ic1,Tw,Tburn,n_sampling, n_neuron,1);
		if low_rate_flag
			error('low firing rate in surrogate simulation!')
		end
    end

	[fa_percentshared, fa_normevals, fa_dshared] = compute_pop_stats(sampling_inds, re, n_neuron, Tw, dim_method);

    single_obj=mean([true_statistics.default_weights(1)*(rate1-true_statistics.rate_mean)^2/true_statistics.rate_var,...
                    true_statistics.default_weights(2)*(FanoFactor1-true_statistics.fano_mean)^2/true_statistics.fano_var]);		
    pair_obj=true_statistics.default_weights(3)*(mean_corr1-true_statistics.mean_corr_mean)^2/true_statistics.mean_corr_var;		
    pop_obj=mean([true_statistics.default_weights(4)*(fa_percentshared-true_statistics.fa_percent_mean)^2/true_statistics.fa_percent_var,...
                        true_statistics.default_weights(5)*(fa_dshared-true_statistics.fa_dim_mean)^2/true_statistics.fa_dim_var,...
                        true_statistics.default_weights(6)*norm(fa_normevals-true_statistics.fa_normeval_mean,2)^2/true_statistics.fa_normeval_var]);		


	
	switch statistics_group
	case '1'
		objective=single_obj/sum(true_statistics.default_weights);
	case '2'
		objective=pair_obj/sum(true_statistics.default_weights);
	case '3'
		objective=pop_obj/sum(true_statistics.default_weights);
	case '12'
		objective=(2*single_obj+pair_obj)/sum(true_statistics.default_weights);
	case '13'
		objective=(2*single_obj+3*pop_obj)/sum(true_statistics.default_weights);
	case '23'
		objective=(pair_obj+3*pop_obj)/sum(true_statistics.default_weights);
	case '123'
		objective=(2*single_obj+pair_obj+3*pop_obj)/sum(true_statistics.default_weights);
	end
catch 
	warning_flag=1;
	warning('something wrong happened in main simulation')
end

if warning_flag
	objective=NaN;
end

if save_stats
	
		[full_stats_to_record] = generate_full_stats_to_record(statistics_group,objective,rate1,var1,FanoFactor1,mean_corr1,fa_percentshared,fa_dshared,fa_normevals);
		full_stats_time=toc;
		execution_time=table(surrogate_time,full_stats_time);
		surrogate_stats=table(cost_surrogate,rate0,var0,FanoFactor0);
		stats_weights=true_statistics.default_weights;
		save_statistics(stats_weights,x,full_stats_to_record,surrogate_stats,execution_time,stats_filename);

end


end