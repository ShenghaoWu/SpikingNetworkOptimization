function [full_stats_to_record] = generate_full_stats_to_record(statistics_group,objective,rate1,var1,FanoFactor1,mean_corr1,fa_percentshared100,fa_dshared100,fa_normevals100)
%Organize statistics into a table for logging
switch statistics_group
	case '1'
		full_stats_to_record=table(objective,rate1,var1,FanoFactor1);
	case '2'
		full_stats_to_record=table(objective,mean_corr1);
	case '3'
		full_stats_to_record=table(objective,fa_percentshared100,fa_dshared100,fa_normevals100);
	case '12'
		full_stats_to_record=table(objective,rate1,var1,FanoFactor1,mean_corr1);

	case '13'
		full_stats_to_record=table(objective,rate1,var1,FanoFactor1,fa_percentshared100,fa_dshared100,fa_normevals100);

	case '23'
		full_stats_to_record=table(objective,mean_corr1,fa_percentshared100,fa_dshared100,fa_normevals100);

	case '123'
		full_stats_to_record=table(objective,rate1,var1,FanoFactor1,mean_corr1,fa_percentshared100,fa_dshared100,fa_normevals100);
end