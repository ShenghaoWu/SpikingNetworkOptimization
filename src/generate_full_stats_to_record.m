function [full_stats_to_record] = generate_full_stats_to_record(statistics_group,objective,rate1,var1,FanoFactor1,mean_corr1,fa_percentshared,fa_dshared,fa_normevals)
%% Organize statistics into a table for logging
%   -Input
%      statistics_group: string, types of statistics in the cost function
%      ('1': single-neuron, '2': pairwise, '3': population)
%      objective: float, cost value
%      rate1,var1,FanoFactor1,mean_corr1,fa_percentshared,fa_dshared,fa_normevals:
%      values of the activity statistics

%   -Output
%      full_stats_to_record: table; activity statistics of the evaluated
%      statistics in a table.


%Organize statistics into a table for logging

switch statistics_group
	case '1'
		full_stats_to_record=table(objective,rate1,var1,FanoFactor1);
	case '2'
		full_stats_to_record=table(objective,mean_corr1);
	case '3'
		full_stats_to_record=table(objective,fa_percentshared,fa_dshared,fa_normevals);
	case '12'
		full_stats_to_record=table(objective,rate1,var1,FanoFactor1,mean_corr1);

	case '13'
		full_stats_to_record=table(objective,rate1,var1,FanoFactor1,fa_percentshared,fa_dshared,fa_normevals);

	case '23'
		full_stats_to_record=table(objective,mean_corr1,fa_percentshared,fa_dshared,fa_normevals);

	case '123'
		full_stats_to_record=table(objective,rate1,var1,FanoFactor1,mean_corr1,fa_percentshared,fa_dshared,fa_normevals);
end