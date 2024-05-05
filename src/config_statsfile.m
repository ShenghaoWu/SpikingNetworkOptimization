function filename=config_statsfile(filename)
%% Initialize an empty statistics logging file
%   -Input
%      filename: string, name of the logging file
%    
%   -Output
%      filename: string, name of the logging file


try 
	load(filename);
	warning('filename exists!!!!!!!!!!!!!!!!!!!!!!!')

catch
	paras=[];
	full_stats=[];
	surrogate_stats=[];
	execution_time=[];
	save(filename,'paras','full_stats','surrogate_stats','execution_time');
end