function [x_star] = random_sample(x_range)
%% Random sample from the ranges
%   -Input
%      x_range: [number of parameters, 2], ranges (lower/upper bounds) for
%      sampling the parameter sets
%     
%   -Output
%      x_star: [number of parameters]; the randomly-sampled parameter set


rng shuffle;

lb = x_range(:,1)';
ub = x_range(:,2)';
x_star = (ub-lb).*rand(1,size(x_range,1)) + lb;


end