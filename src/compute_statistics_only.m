function [rate1,var1, FanoFactor1, mean_corr1]=compute_statistics_only(re1_s) 
%% Compute the rate, var, Fano, and rsc of the spike count matrix.

%   -Input
%      re1_s: [number of neurons, number of bins], spike count matrix
%      
%   -Output
%      rate1: float;  mean firing rate
%      var1: float, mean variance
%      FanoFactor1: float, mean Fanofactor
%      mean_corr1: float, mean rsc



rng shuffle;
COV=cov([re1_s]');
Var=diag(COV);
var1=mean(Var);
rate1=mean(mean(re1_s,2));
FanoFactor1 = mean(Var./mean(re1_s,2));
R = COV./sqrt(Var*Var'); 
upper_R=R(triu(true(size(R)), 1));
mean_corr1=nanmean(rtoZ(upper_R));
end


function Z = rtoZ(r)
% Z = rtoZ(r)
%
% RTOZ translates fisher r correlations into Z scores
%  aka the "Fisher r-to-Z' transformation"

% Matthew A. Smith
% Revised: 2001.04.21
% Update: nan is permitted here
% jittering when r = 1
if (~isempty(find(r<=-1-1e-10,1)) || ~isempty(find(r>=1+1e-10,1))) & ~isnan(r)
    error('r values must be bounded by -1 and 1');
end

Z = 0.5 * log((1 + r+1e-10) ./ (1 - (r-1e-10)));
end