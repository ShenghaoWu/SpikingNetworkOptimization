function [obj_configs, optimization_opt] =  generate_config(varargin)
% Generate the configuration struct for SNOPS
% To specify a configuration field, input the name,value pair to the function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Configs for SNN simulation%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


is_save=0; % save data
CompCorr=0; % compute correlations
Layer1only=1; % 1 for two-layer network, 0 for three-layer network
loadS1=0; % if load pre-computed ffwd spiketrains
plotPopR=0; % plot population rate
fixW=0;  % use the same weight matrices for multiple simulations
verbose=0; % if print info during SNN simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Configs for evaluating the objective(cost) function%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = ''; % name of the log file to be saved throughout customization, should end with "_bo_output"
root = './data/'; % root folder for the target datasets
save_root='./results/'; % root folder for the save files
real_data_name = ''; %name of the file that contains the target activity statistics 
is_surrogate=1; %if running a short simulation to determine feasibility
n_sampling=10; %number of samplings of the neurons for computing the activity statistics
Tw=200; % time window for computing the activity statistics
Tburn=500; % time window at the beginning to exclude (remove transcient state)
T0=10000; %simulation length (ms) for the short simulation
T_short=20000;  %simulation length (ms) for the full simulation when only computing single/pairwise stats
T=700*200+500; %simulation length (ms) for the full simulation when computing single/pairwise/pop stats
dt=0.05; % simulation time step
N_groundtruth=50000; % size of the original network from Huang et al, 2019, used for scaling to a smaller size
tolerance=Inf; % deprecated; threshold for the cost of the small newtork simulation
metric_norm='L2'; % difference function for the cost
statistics_group='123'; % which of single (1) pairwise (2) pop (3) stats to include
survival_rate=Inf; % deprecated
is_simulation=1; % deprecated
is_small=1; % use the smaller network size for network simulation
Ne1=50; % number of E neurons in the recurrent layer per side, actual number is Ne1^2
Ni1=25; % number of I neurons in the recurrent layer per side, actual number is Ne1^2
save_stats=1; % if save stats during customization
is_spatial=1; % if use the SBN, if not, use the CBN
dim_method = 'CV_skip'; % cross-validation method for FA, skipping some CV dims to speed up
stats_weights = [1,1,1,1,1,1]; % weights for each statistic.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Configs for Bayesian Optimization%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_range=[1,25;1,25;0,0.25;0,0.25;0,0.25;-150,0;0,150;-150,0;0,150;0,150;0,150]; % parameter search range (lb,ub)
epsilon=0; % controls exploration ratio in the acquisition function, can be tuned but we simply use 0 here.
n_sample=100000; % number of total samples to draw for optimization the acquistion function
n_local=10; % number of local samples to draw for optimization the acquistion function
n_data=50; %number of the initial random samples for BO
max_iter=10000; % max number of iterations
max_time=10000*3600; % max customization time
std_tol=0.15; % std tolerance for intensification
min_cost_eval=5; % min number of repeated evaluations for the same parameter set in intensification
max_cost_eval=10; % max number of repeated evaluations for the same parameter set in intensification
n_check=1e6; % iteration interval for checking the log files of concurrent threads for results recycling
is_log=1; % if using log(cost) for fitting GP
base_name = []; % if not empty, use the base network simulation samples as the initial seed points
top_n=50; % number of initial seed samples for the base file
sort_mode='random'; % sorting method for the base file
is_sort=0; % if sort the base file
cost_stopping_thres=0.1; %threshold of the cost below which to stop SNOPS.

%assign values from varagin
for cnt = 1: length(varargin)/2
  eval([varargin{cnt*2-1} '=varargin{cnt*2};']);
end


opt={};
obj_configs={};
optimization_opt={};

%%%%%%%%%%%%Assigning options for SNN%%%%%%%%%%%%%%%%%%%%%%%
opt.save=is_save; 
opt.CompCorr=CompCorr;
opt.Layer1only=Layer1only; 
opt.loadS1=loadS1;
opt.plotPopR=plotPopR;
opt.fixW=fixW; 
opt.verbose=verbose;

%%%%%%%%%%%%Assigning options for cost%%%%%%%%%%%%%%%%%%%%%%%
obj_configs.opt=opt;
obj_configs.filename=filename;
obj_configs.root=root;
obj_configs.save_root=save_root;
obj_configs.is_surrogate=is_surrogate;

obj_configs.n_sampling=n_sampling;
obj_configs.Tw=Tw; 
obj_configs.Tburn=Tburn; 
obj_configs.T0=T0; 
obj_configs.T_short=T_short; 
obj_configs.T=T;
obj_configs.dt=dt;
obj_configs.N_groundtruth=N_groundtruth;
obj_configs.tolerance=tolerance;
obj_configs.metric_norm=metric_norm;
obj_configs.statistics_group=statistics_group;
obj_configs.survival_rate=survival_rate;
obj_configs.is_simulation=is_simulation;
obj_configs.is_small=is_small;
obj_configs.Ne1=Ne1;
obj_configs.Ni1=Ni1;
obj_configs.save_stats=save_stats;
obj_configs.is_spatial=is_spatial;
obj_configs.dim_method = dim_method;
obj_configs.stats_weights = stats_weights;

%%%%%%%%%%%%Assigning options for optimization%%%%%%%%%%%%%%%%%%%%%%%
optimization_opt.x_range=x_range;
optimization_opt.epsilon=epsilon;
optimization_opt.n_sample=n_sample;
optimization_opt.n_local=n_local;
optimization_opt.n_data=n_data;
optimization_opt.max_iter=max_iter;
optimization_opt.max_time=max_time;
optimization_opt.std_tol=std_tol;
optimization_opt.min_cost_eval=min_cost_eval;
optimization_opt.max_cost_eval=max_cost_eval;
optimization_opt.n_check=n_check;
optimization_opt.is_log=is_log;
optimization_opt.base_name = strcat(save_root,base_name,'.mat');
optimization_opt.top_n=top_n;
optimization_opt.sort_mode=sort_mode;
optimization_opt.is_sort=is_sort;
optimization_opt.cost_stopping_thres=cost_stopping_thres;

  

%configure log filenames
filename=obj_configs.filename; %filename contains jobid!
results_name=strcat(obj_configs.save_root,filename,'.mat');
obj_configs.real_data_name = strcat(root,real_data_name,'.mat');
load(obj_configs.real_data_name);

obj_configs.true_statistics=true_statistics;


stats_name=strcat(obj_configs.save_root,filename,'_stats.mat');
obj_configs.stats_filename=config_statsfile(stats_name);
optimization_opt.save_name=results_name;

%configure the target (true) statistics
obj_configs.true_statistics.default_weights=obj_configs.stats_weights;

%configure base for intensifications
if ~isempty(base_name)
	[base_name] = base_real(obj_configs.real_data_name,optimization_opt.base_name,obj_configs.filename,optimization_opt.top_n,optimization_opt.sort_mode,obj_configs.stats_weights,optimization_opt.is_sort);
    optimization_opt.base_name = base_name;
else 
	optimization_opt.base_name = [];
end




end
