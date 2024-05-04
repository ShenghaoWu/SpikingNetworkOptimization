function [output] = EI_interaction(x,gp_object,gp_object_feas,f_plus,epsilon)

%x = cat(2,x(:,:),x(:,1)./x(:,2),x(:,3)./x(:,4),x(:,6).*x(:,7));
x = cat(2,x(:,:),x(:,3)./x(:,4),x(:,6).*x(:,7)-x(:,8).*x(:,9));


[miu_x,sigma_x]= predict(gp_object,x);
[mu_feas,std_feas]= predict(gp_object_feas,x);
prob_feas = normcdf((mu_feas-.5)./std_feas);
Z=(f_plus-epsilon-miu_x)./sigma_x;
output=prob_feas.*((f_plus-epsilon-miu_x).*normcdf(Z)+sigma_x.*normpdf(Z));
