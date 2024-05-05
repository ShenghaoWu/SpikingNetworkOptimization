function [unstable_flag]=is_unstable(re) 
%% Feasibility constraint checker for unstable solutions using changed point detection.
%% Check if the mean of first 5 entires is within 3 std of the mean of the last 5 entries
%%       to rule out the case of long transient oscillation
%   -Input
%      re: [number of neurons, number of bins], spike count matrix

%   -Output
%      unstable_flag: {0,1}; if the spiking activity is unstable (long
%      transcient or switching dynamics).
ipt = findchangepts(mean(re,1));
unstable_flag = abs(mean(mean(re(:,1:ipt),1))-mean(mean(re(:,ipt+1:end),1)))>=3*std(mean(re(:,ipt+1:end),1));

