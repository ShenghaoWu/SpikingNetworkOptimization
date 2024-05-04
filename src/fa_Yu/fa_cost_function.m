function [costs] = fa_cost_function(paras, real_eigs, real_pshared,varargin)
%
% [costs] = fa_cost_function(paras, ...) 
%
% Calculate the cost given the parameters
%
% INPUTS: 
%
% paras    - parameters
%
% OUTPUTS:
%
% costs  - costs


%X = simdata_fa(paras, 5000);
[eig_vectors, eigs]=eig(paras.L*paras.L');
[eigs,idx]=sort(diag(eigs),'descend');


[percentshared, d_shared,normevals] = compute_shared(paras, 0.95);
eigs_norm = norm(eigs-real_eigs,1)/norm(real_eigs,1);
pshared_norm=percentshared/real_pshared;
costs=0.01*eigs_norm+0.99*pshared_norm;                                  
%costs=eigs_norm;