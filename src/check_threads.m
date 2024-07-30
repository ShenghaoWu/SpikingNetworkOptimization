function [x_star,y_star,is_better]=check_threads(x_train,y_train,root,stats_name,true_statistics) 
%% Function for results recycling. This function  will load all the stats log files of the same name patterns (middle_str) from 1 to nmax and compute the cost based on saved activity statistics and the target stats (true_statistics).
%%	If we encounter a network simulation yielding lower cost than the current best, return it.
%   -Input
%		x_train: [number of samples, number of parameters]
%		y_train: [number of samples, 1]
%		root: string; folder containing logging files 
% 		stats_name: string; name of the stats logging files
%		true_statistics: table; target statistics
%   -Output
%       x_star:  [1, number of parameters] the potential optimal parmaeter set, if exists. 
%       y_star:  [1] cost of the potential optimal parmaeter set, if exists.
%       is_better:  {0,1}, if the potential optimal parmaeter exists.



x_star=nan;
y_star=nan;
is_better=0;
nmax= 100;
middle_str =  extractBetween(stats_name,root,'_output_');
stats = [];
ps = [];





for jobid=1:nmax
    try
        load(strcat(root,middle_str,'_output_',string(jobid),'_stats.mat'));
        stats=[stats;full_stats];
        ps=[ps;paras];
    catch 
        
    end
    
end



full_stats=stats;

rate1_cost=true_statistics.default_weights(1)*(full_stats.rate1-true_statistics.rate_mean).^2/true_statistics.rate_var;
FanoFactor1_cost = true_statistics.default_weights(2)*(full_stats.FanoFactor1-true_statistics.fano_mean).^2/true_statistics.fano_var;
mean_corr1_cost = true_statistics.default_weights(3)*(full_stats.mean_corr1-true_statistics.mean_corr_mean).^2/true_statistics.mean_corr_var;
fa_percentshared100_cost = true_statistics.default_weights(4)*(full_stats.fa_percentshared100-true_statistics.fa_percent_mean).^2/true_statistics.fa_percent_var;
fa_dshared100_cost = true_statistics.default_weights(5)*(full_stats.fa_dshared100-true_statistics.fa_dim_mean).^2/true_statistics.fa_dim_var;
fa_normevals100_cost=true_statistics.default_weights(6)*vecnorm(full_stats.fa_normevals100-true_statistics.fa_normeval_mean,2,2).^2/true_statistics.fa_normeval_var;
objs=sum([rate1_cost,FanoFactor1_cost,mean_corr1_cost,fa_percentshared100_cost,fa_dshared100_cost,fa_normevals100_cost],2)/sum(true_statistics.default_weights);


logicals = ~ismember(ps{:,1}, x_train(:,1));

if sum(logicals)>0
	objs = objs(logicals);
	ps = ps{logicals,:};
	[I,J]=min(objs);
	if (isnan(nanmin(y_train)) && ~isnan(I))|| I<nanmin(y_train)
        is_better=1;
        x_star=ps(J,:);
        y_star=I;

	end 
end


end