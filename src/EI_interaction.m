function [output] = EI_interaction(x,gp_object,gp_object_feas,f_plus,epsilon)

%% Function for computing the acquisition using the expectation maximization formula when accounting for parameter interactions for the timescales and connectivity strengths
%% Considering both the full cost GP and the feasibility GP
%   -Input
%      x: [number of parameters], numeric parameter set
%      gp_object: struct, fitted GP for the full cost function
%      gp_object: struct, fitted GP for the feasibility constraints
%      f_plus: float, current best cost value
%      epsilon: float, exploration-exploitation trade-off parameter

%   -Output
%      output: float; the EI value

x = cat(2,x(:,:),x(:,3)./x(:,4),x(:,6).*x(:,7)-x(:,8).*x(:,9));

[miu_x,sigma_x]= predict(gp_object,x);
[mu_feas,std_feas]= predict(gp_object_feas,x);
prob_feas = normcdf((mu_feas-.5)./std_feas);
Z=(f_plus-epsilon-miu_x)./sigma_x;
output=prob_feas.*((f_plus-epsilon-miu_x).*normcdf(Z)+sigma_x.*normpdf(Z));
