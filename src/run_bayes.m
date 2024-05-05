function run_bayes(obj_configs, optimization_opt, optimization_algo)
%% Entry point for the SNOPS customization.
%   -Input
%		obj_configs: struct, configurations for the network simulation
%		optimization_opt: struct, configurations for the Bayesian optimization
%       optimization_algo: int, 1:use BO (SNOPS default); 
%                               2: use accelerated random search;
%                               3: use the typical random search (without feasibility constraints and intensification);
%                               4: use BO with interaction terms (during GP
%                               fitting).



%Random seed for running on cluster
seed_offset = randi(floor(intmax/10));
rng(rand(1)*1000 + seed_offset);
warning('optmization begins!')
if isfile(optimization_opt.save_name) %resume if results files already exist
	warning('file exists, resuming!')
	load(optimization_opt.save_name)
	optimization_opt.x_train=x_train;
	optimization_opt.y_train=y_train;
	optimization_opt.incumbent_std=incumbent_std;
	optimization_opt.optimization_time=optimization_time;
	optimization_opt.y_feasibility=y_feasibility;
elseif ~isempty(optimization_opt.base_name) %intensify if base file is there; use intensification_vanilla for grid search
	warning('intensifying base results!')
	if is_random==2
		[optimization_opt,obj_configs] = intensification_vanilla(optimization_opt.base_name,obj_configs,optimization_opt);
	else
		[optimization_opt,obj_configs] = intensification(optimization_opt.base_name,obj_configs,optimization_opt);
	end
	
end

%Preparing cost function wrapper
fun=@(x,is_surrogate)cost_func(x,is_surrogate,obj_configs);
func=@(x,is_surrogate)cost_interface(x,is_surrogate,fun);	

switch optimization_algo

case 1
	fprintf('Using Bayesian optimization as the backbone, customization begins.')
	[x_train,y_train,y_feasibility,optimization_time] = bayesian_optimization(func,optimization_opt,obj_configs);

case 2
	fprintf('Using accelerated random search as the backbone, customization begins.')
	[x_train,y_train,y_feasibility,optimization_time] = random_search(func,optimization_opt,obj_configs); 

case 3
	fprintf('Using random search as the backbone, customization begins.')
	[x_train,y_train,y_feasibility,optimization_time] = random_search_vanilla(func,optimization_opt,obj_configs); 

case 4
	fprintf('Using Bayesian optimization with parameter interaction as the backbone, customization begins.')
	[x_train,y_train,y_feasibility,optimization_time] = bayesian_optimization_interaction(func,optimization_opt,obj_configs);

otherwise
	error('Unrecognized optimization algo id!')
end

end

