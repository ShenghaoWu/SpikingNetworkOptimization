function [gprMdl,gprMdl_feas,f_plus,epsilon_scaled]=fit_gp(x_train,y_train,y_feasibility,epsilon,gprMdl,current_feas,f_plus,epsilon_scaled,is_log) 
%Fit the GP models on the evaluated parameter sets, for both the full cost and the feasibility constraints.

if is_log
    y_train=log(y_train);
end

if current_feas
    try
        gprMdl = fitrgp(x_train,y_train,...
                'KernelFunction','ardmatern52',...
                'BasisFunction', 'constant', ...
                'Standardize', false,...
                'ConstantSigma',false,...
                'FitMethod','exact','PredictMethod','exact');
    catch 
        % If no convergence, set a lower bound for sigma in the kernel
        gprMdl = fitrgp(x_train,y_train,...
                'KernelFunction','ardmatern52',...
                'BasisFunction', 'constant', ...
                'Standardize', false,...
                'ConstantSigma',false,...
                'SigmaLowerBound',1e-5,...
                'FitMethod','exact','PredictMethod','exact');
    end
    [mu_train,~]= predict(gprMdl,x_train(find(y_feasibility==1),:));
    [f_plus,~]=min(mu_train);
    epsilon_scaled=nanstd(y_train)*epsilon;
end

try
    gprMdl_feas = fitrgp(x_train,y_feasibility,...
            'KernelFunction','ardmatern52',...
            'BasisFunction', 'constant', ...
            'Standardize', false,...
            'ConstantSigma',false,...
            'FitMethod','exact','PredictMethod','exact');
catch 
    % If no convergence, set a lower bound for sigma in the kernel
    gprMdl_feas = fitrgp(x_train,y_feasibility,...
            'KernelFunction','ardmatern52',...
            'BasisFunction', 'constant', ...
            'Standardize', false,...
            'ConstantSigma',false,...
            'SigmaLowerBound',1e-5,...
            'FitMethod','exact','PredictMethod','exact');
end

end