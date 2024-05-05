function [rate1,var1,FanoFactor1,mean_corr1,fa_percentshared,fa_dshared,fa_normevals] = nan_stats(m)
%% Return NaN-filled activiy statistics
%   -Input
%      m: int, dimensionality for the eigenspectrum
%   -Output
%      rate1,var1,FanoFactor1,mean_corr1,fa_percentshared,fa_dshared,fa_normevals,
%    : nan-filled activity statistics.

rate1=NaN;
var1=NaN;
FanoFactor1=NaN;
mean_corr1=NaN;
fa_percentshared=NaN;
fa_dshared=NaN;
fa_normevals=NaN(1,m);
