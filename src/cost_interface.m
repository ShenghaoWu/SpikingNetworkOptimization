function [result] = cost_interface(x,is_surrogate,fun)
%% Wrapper function which takes the parameter vector, converts it to a named table, and feed it to the cost function.
%   -Input
%       x: [number of parameters], numeric parameter set for the simulation
%		is_surrogate: 0 (don't run the short simulation for feasibility), 1 (run the short simulation)
% 		fun: function, the actual cost function wrapper (e.g. cost_func.m)
%   -Output
%       result: float, cost value corresponding to x.

if length(x) == 8 % CBN
	taudsynI=x(1);
	taudsynE=x(2);
	JrEI= x(3);
	JrIE= x(4);
	JrII= x(5);
	JrEE= x(6);
	JrEX= x(7);
	JrIX= x(8);
	input_paras=table(taudsynI,taudsynE,JrEI, JrIE, JrII, JrEE,JrEX,JrIX);
	result=fun(input_paras,is_surrogate);

elseif length(x) == 5 % CBN, fewer parameters 
	taudsynI=x(1);
	mean_sigmaRRIs=x(2);
	JrEI= x(3);
	JrIE= x(4);
	JrEX= x(5);
	input_paras=table(taudsynI,mean_sigmaRRIs,JrEI, JrIE,JrEX);
	result=fun(input_paras,is_surrogate);

else
	% SBN, fewer parameters 

	taudsynI=x(1);
	taudsynE=x(2);
	mean_sigmaRRIs=x(3);
	mean_sigmaRREs=x(4);
	mean_sigmaRXs= x(5);
	JrEI= x(6);
	JrIE= x(7);
	JrII= x(8);
	JrEE= x(9);
	JrEX= x(10);
	JrIX=x(11);
	input_paras=table(taudsynI,taudsynE,mean_sigmaRRIs,mean_sigmaRREs, mean_sigmaRXs, JrEI, JrIE, JrII, JrEE,JrEX,JrIX);
	result=fun(input_paras,is_surrogate);


end
end