function [output] = EI(x,gp_object,gp_object_feas,f_plus,epsilon)
%% Function for computing the acquisition using the expectation maximization formula
%% Considering both the full cost GP and the feasibility GP
%   -Input
%      x: [number of parameters], numeric parameter set
%      gp_object: struct, fitted GP for the full cost function
%      gp_object_feas: struct, fitted GP for the feasibility constraints
%      f_plus: float, current best cost value
%      epsilon: float, exploration-exploitation trade-off parameter

%   -Output
%      output: float; the EI value

[miu_x,sigma_x]= predict(gp_object,x);
[mu_feas,std_feas]= predict(gp_object_feas,x);
prob_feas = normcdf((mu_feas-.5)./std_feas);
Z=(f_plus-epsilon-miu_x)./sigma_x;
output=prob_feas.*((f_plus-epsilon-miu_x).*normcdf(Z)+sigma_x.*normpdf(Z));
