function [optimization_opt,obj_configs] = intensification_vanilla(base_name,obj_configs,optimization_opt)
%% Intensification for the initial base file, only if the base file exists. For the random search method (see intensification.m for BO).
%   -Input
% 		base_name: string, filename for the previously generated network
% 		simulations.
%		obj_configs: struct, configurations for the network simulation
%		optimization_opt: struct, configurations for the Bayesian optimization

%   -Output
%		optimization_opt: struct, updated configurations for the Bayesian optimization
%		obj_configs: struct, updated configurations for the network simulation


load(base_name)
paras_=paras;
full_stats_=full_stats;
surrogate_stats_=surrogate_stats;
execution_time_=execution_time;
load(obj_configs.stats_filename);
paras=[paras;paras_];
full_stats=[full_stats;full_stats_];
surrogate_stats=[surrogate_stats;surrogate_stats_];
execution_time=[execution_time;execution_time_];
save(obj_configs.stats_filename,'paras','full_stats','surrogate_stats','execution_time');
load(base_name)
true_statistics= obj_configs.true_statistics;




rate1_cost=true_statistics.default_weights(1)*(full_stats.rate1-true_statistics.rate_mean).^2/true_statistics.rate_var;
FanoFactor1_cost = true_statistics.default_weights(2)*(full_stats.FanoFactor1-true_statistics.fano_mean).^2/true_statistics.fano_var;
mean_corr1_cost = true_statistics.default_weights(3)*(full_stats.mean_corr1-true_statistics.mean_corr_mean).^2/true_statistics.mean_corr_var;
fa_percentshared_cost = true_statistics.default_weights(4)*(full_stats.fa_percentshared-true_statistics.fa_percent_mean).^2/true_statistics.fa_percent_var;
fa_dshared_cost = true_statistics.default_weights(5)*(full_stats.fa_dshared-true_statistics.fa_dim_mean).^2/true_statistics.fa_dim_var;
fa_normevals_cost=true_statistics.default_weights(6)*vecnorm(full_stats.fa_normevals-true_statistics.fa_normeval_mean,2,2).^2/true_statistics.fa_normeval_var;
full_stats{:,1}=mean([rate1_cost,FanoFactor1_cost,mean_corr1_cost,fa_percentshared_cost,fa_dshared_cost,fa_normevals_cost],2);

if obj_configs.is_spatial
	x_train=paras{:,:};
else
	x_train=[paras.taudsynI,paras.taudsynE,paras.JrEI,paras.JrIE,paras.JrII,paras.JrEE,paras.JrEX,paras.JrIX];
end
y_train_eval=full_stats.objective;


fun=@(x,is_surrogate)cost_func(x,is_surrogate,obj_configs);
func=@(x,is_surrogate)cost_interface(x,is_surrogate,fun);	



optimization_time=execution_time.surrogate_time+execution_time.full_stats_time;

n_data=length(optimization_time);
y_feasibility=NaN(n_data,1);
y_train=NaN(n_data,1);
try
	incumbent_std=optimization_opt.incumbent_std;
	current_min = optimization_opt.current_min;
catch 
	incumbent_std=0;
	current_min=nanmin(y_train);
end

max_cost_eval=optimization_opt.max_cost_eval;
min_cost_eval=optimization_opt.min_cost_eval;
std_tol=optimization_opt.std_tol;
for i =1:n_data
    t_init=tic;
    y_hat=y_train_eval(i);
    tmp_yhat=[y_hat];
    for kkk = 1:4
        tmp_yhat=[tmp_yhat,func(x_train(i,:),obj_configs.is_surrogate)];
    end
    y_hat=mean(tmp_yhat);
    y_feasibility(i)=0;
    y_train(i)=y_hat;
    current_feas=0;
    optimization_time(i)=optimization_time(i)+toc(t_init);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%determine threshold%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if obj_configs.survival_rate == Inf 
	tolerance = Inf ;
else
	rate0_cost=true_statistics.default_weights(1)*(surrogate_stats.rate0-true_statistics.rate_mean).^2/true_statistics.rate_var;
	FanoFactor0_cost = true_statistics.default_weights(2)*(surrogate_stats.FanoFactor0-true_statistics.fano_mean).^2/true_statistics.fano_var;
	surrogate_obj=mean([rate0_cost,FanoFactor0_cost],2);
	survivors = sort(surrogate_obj(~isnan(surrogate_obj)));
	for scnt=1:length(survivors)
	    survival_rate=sum(surrogate_obj<=survivors(scnt))/length(surrogate_obj);
	    if survival_rate>=obj_configs.survival_rate
	        tolerance = survivors(scnt);
	        break;
	    end
	    tolerance = survivors(scnt);
	end
end

obj_configs.tolerance=tolerance;
name_to_save=obj_configs.stats_filename;
name_to_save=strcat(erase(name_to_save,'.mat'),'_intens.mat');
save(name_to_save,'tolerance','current_min','incumbent_std');




optimization_opt.x_train=x_train;
optimization_opt.y_train=y_train;
n_data=size(x_train,1);
optimization_opt.n_data = n_data;
optimization_opt.incumbent_std=incumbent_std;
optimization_opt.optimization_time=optimization_time;
optimization_opt.y_feasibility=y_feasibility;

