%% Full demonstration of SNOPS (running time ~ 168 hours)
% This script demonstrates the customization of a SBN to 10 example simulation 
% datasets generated  by the SBN with varying parameter sets.
% Make sure you cd to /SpikingNetworkOptimization, 
% then execute the following scripts to generate the config files and start SNOPS. 
% For details on the customization configurations, please refer to the documentation of generate_config.m
% The code below assumes you are running on a machine with at least 10 CPU cores, 
% as it will spawn  10 batch jobs in the background and use the results recycling mechanism 
% which routinely checks the saved results from the concurrent threads. 
% Expect 7 days to finish the full customization (slurm/screen/tmux recommended), but you may
% be able to obtain reasonable results for some of the target datasets after  48-72 hours. 
% If you do not want to use results recycling, simply change n_check to 1e6, which corresponds to 
% the routine interval for the recycling. This may possibly yield slower convergence or suboptimal results.
% You can load the log file under the results folder with the script below after or during the customization.


clear; clc; close all;
addpath('src')


for i = 1:10
    [obj_configs,optimization_opt]=generate_config('is_spatial', 1,... % =1 using the SBN, =0 using the CBN.
                                                   'max_time', 3600*7*24,... % max running time (in seconds).
                                                   'n_check', 10,... % results recycling checking interval. if changed to 1e6, turn off recycling.
                                                   'min_cost_eval', 2,... %number of min evaluations for intensification, the larger the more stable final results are.
                                                   'max_cost_eval', 3,... %number of max evaluations for intensification, the larger the more stable final results are.
                                                   'x_range',[1,25; 1,25; 0,0.25; 0,0.25; 0,0.25; -150,0; 0,150; -150,0; 0,150; 0,150; 0,150],... % search region, in [n_params, 2]. For the SBN, this  correspond to: taudsynI, taudsynE, mean_sigmaRRIs, mean_sigmaRREs, mean_sigmaRXs, JrEI, JrIE, JrII, JrEE, JrEX, JrIX
                                                   'real_data_name', strcat('demo_sbn_simu_',string(i)),... % name pattern for the target stats, must be placed under /data.
                                                   'filename', strcat('demo_recy_bo_output',string(i))); % name pattern for the logging files, will be placed under /results. Note: filename must end with "_bo_output" for the recyling mechanism to identify all concurrent files

    optimization_algo=1; % Use Bayesian optimization as the backbone
    eval(strcat(obj_configs.filename,'=batch(@()run_bayes(obj_configs, optimization_opt, optimization_algo), 0, {});'));

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%loading logging files (uncomment & run when you want to check the results)%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
file_pattr = 'demo_recy_bo_output';
nthread = 10;
min_costs=[]; %minimal cost for each customization task
optimal_paras=[]; %optimal parameter set for each customization task
min_stats_mean = zeros(nthread, 55);  %activity stats of the optimal parameter set (es has 50 entries, leading to a total of 55)
target_stats_mean = zeros(nthread, 55); %activity stats of the target data
cost_trace = {}; %trace of the cost over iteration for each customization task
time_trace = {}; %trace of the time over iteration for each customization task

for jobid=1:nthread
  y_trains = [];
  x_trains = [];
  stats = [];
  parass = [];
  try
    results_name = strcat('./results/', file_pattr, string(jobid), '.mat');
    stats_name=strcat('./results/', file_pattr, string(jobid), '_stats.mat');
    load(results_name)
    load(stats_name)
    y_trains=[y_trains;y_train];
    x_trains=[x_trains;x_train];
    stats = [stats;full_stats];
    parass = [parass;paras];
    cost_trace{end+1} = y_train;
    time_trace{end+1} = optimization_time;
    [I, J]=min(y_trains);
    min_costs = [min_costs, I];
    optimal_paras = [optimal_paras; x_trains(J, :)];
    pa1 = x_trains(J, 1);
    J = parass{:, 1} == pa1;
    pas = parass{J, :};
    pas = pas(1, :);
    min_stats_mean(jobid, :)=mean(stats{J,[2,4:end]},1);

    load(strcat('./data/demo_sbn_simu_',string(jobid),'.mat'))
    target_stats_mean(jobid,:) = [true_statistics.rate_mean, true_statistics.fano_mean, true_statistics.mean_corr_mean, true_statistics.fa_percent_mean, true_statistics.fa_dim_mean,true_statistics.fa_normeval_mean];


    fprintf('successfully loaded file %d \n', jobid)
    fprintf('target stats: fr: %.2f, ff: %.2f, rsc: %.3f, psh: %.1f, dsh: %.1f, es (1st): %.1f \n',...
            target_stats_mean(jobid,1), target_stats_mean(jobid,2), target_stats_mean(jobid,3),...
            target_stats_mean(jobid,4)*100, target_stats_mean(jobid,5), target_stats_mean(jobid,6) )
    fprintf('customized stats: fr: %.2f, ff: %.2f, rsc: %.3f, psh: %.1f, dsh: %.1f, es (1st): %.1f \n',...
            min_stats_mean(jobid,1), min_stats_mean(jobid,2), min_stats_mean(jobid,3),...
            min_stats_mean(jobid,4)*100, min_stats_mean(jobid,5), min_stats_mean(jobid,6) )
    fprintf('optimal parameter set: taudsynI: %.2f, taudsynE: %.2f, mean_sigmaRRIs: %.2f, mean_sigmaRREs: %.2f, mean_sigmaRXs: %.2f, JrEI: %.2f, JrIE: %.2f, JrII: %.2f, JrEE: %.2f, JrEX: %.2f, JrIX: %.2f \n',...
            optimal_paras(jobid,:))
    disp('~~~~~~~~~~~~~~~~~~~~')

  catch 
    fprintf('error loading file %d\n', jobid)
  end
end
%}

