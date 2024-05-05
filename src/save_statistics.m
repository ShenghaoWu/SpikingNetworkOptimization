function  save_statistics(stats_weights,paras_tmp,full_stats_tmp,surrogate_stats_tmp,execution_time_tmp,stats_filename)
%% Append the statistics, parameters, execution time to the existing statistics log file. 
%   -Input
%     stats_weights: [6], weighting for each activity statistics
%     paras_tmp: table, parameter set to append
%     full_stats_tmp: table, activity statistics of the full simulation to append
%     surrogate_stats_tmp, table, activity statistics of the surrogate simulation to append
%     execution_time_tmp: table, running time for both short and full
%     simulation
%     stats_filename: string, name of the statistics logging file



try
	load(stats_filename);
	paras=[paras; paras_tmp];
	full_stats=[full_stats; full_stats_tmp];
	surrogate_stats=[surrogate_stats; surrogate_stats_tmp];
	execution_time=[execution_time; execution_time_tmp];
	save(stats_filename,'stats_weights','paras','full_stats','surrogate_stats','execution_time');
catch
	% In case of IO conflicts, do it again.
	try
		load(stats_filename);
		paras=[paras; paras_tmp];
		full_stats=[full_stats; full_stats_tmp];
		surrogate_stats=[surrogate_stats; surrogate_stats_tmp];
		execution_time=[execution_time; execution_time_tmp];
		save(stats_filename,'stats_weights','paras','full_stats','surrogate_stats','execution_time');
	
	catch 
	end

end

end
