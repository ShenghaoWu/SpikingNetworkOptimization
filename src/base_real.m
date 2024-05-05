function [base_name] = base_real(real_data_name,base_name,filename,top_n,sort_mode,stats_weights,is_sort)
%% Function for loading previously-evaluated simulations.

%   Input
%       - real_data_name: string, filename of the target statistics.
%       - base_name: string, filename of the previously generated network
%       simulations
%       - filename: string, logging file name.
%       - top_n: int, number of simulations of lowest cost to be selected
%       - sort_mode: string, sorting order.
%       - stats_weights: [6], weightings for each statistic.
%       - is_sort: {0,1}, if sort or use random ordering.

%   Output
%       - base_name: string, filename of the previously generated network.


try 
	load(strcat('./q/',filename,'_base.mat'));
	base_name=strcat('./q/',filename,'_base.mat');
	warning('base already exists')

catch

	load(real_data_name);
	load(base_name);

	rate1_cost=stats_weights(1)*(full_stats.rate1-true_statistics.rate_mean).^2/true_statistics.rate_var;
	FanoFactor1_cost = stats_weights(2)*(full_stats.FanoFactor1-true_statistics.fano_mean).^2/true_statistics.fano_var;
	mean_corr1_cost = stats_weights(3)*(full_stats.mean_corr1-true_statistics.mean_corr_mean).^2/true_statistics.mean_corr_var;
	fa_percentshared_cost = stats_weights(4)*(full_stats.fa_percentshared-true_statistics.fa_percent_mean).^2/true_statistics.fa_percent_var;
	fa_dshared_cost = stats_weights(5)*(full_stats.fa_dshared-true_statistics.fa_dim_mean).^2/true_statistics.fa_dim_var;
	fa_normevals_cost=stats_weights(6)*vecnorm(full_stats.fa_normevals-true_statistics.fa_normeval_mean,2,2).^2/true_statistics.fa_normeval_var;
	full_stats{:,1}=mean([rate1_cost,FanoFactor1_cost,mean_corr1_cost,fa_percentshared_cost,fa_dshared_cost,fa_normevals_cost],2);



	if is_sort
		[~,I] = sort(full_stats{:,1});
	else
		I = randperm(size(full_stats,1));
	end

	I=I(1:top_n);
	paras=paras(I,:);
	full_stats=full_stats(I,:);
	surrogate_stats=surrogate_stats(I,:);
	execution_time=execution_time(I,:);

	switch (sort_mode)
		case 'ascend'
		
		%default is acsend

		case 'descend'

		[~,I] = sort(full_stats{:,1},'descend');
		paras=paras(I,:);
		full_stats=full_stats(I,:);
		surrogate_stats=surrogate_stats(I,:);
		execution_time=execution_time(I,:);
			
		case 'random'
		I=randperm(top_n);
		paras=paras(I,:);
		full_stats=full_stats(I,:);
		surrogate_stats=surrogate_stats(I,:);
		execution_time=execution_time(I,:);


		case 'random_subset' % randomly pick top 100 to evaluate
		I=randperm(top_n);
		paras=paras(I(1:ceil(top_n/10)),:);
		full_stats=full_stats(I(1:ceil(top_n/10)),:);
		surrogate_stats=surrogate_stats(I(1:ceil(top_n/10)),:);
		execution_time=execution_time(I(1:ceil(top_n/10)),:);


	end


	base_name= strcat('./q/',filename,'_base.mat');
	save(base_name,'paras','full_stats','surrogate_stats','execution_time')

end