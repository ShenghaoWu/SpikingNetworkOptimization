function [Ic1] = sample_e_neurons(s1,Ne1,is_simulation)
%% Obtain indices for the recurrent E neurons to compute activity statistics
%   -Input
%     s1: [number spikes, 2], spike trains
%     Ne1: int, number of E neurons in the recurrent layer per side (total
%     number if Ne1^2)
%     is_simulation: {0,1}, if = 1, does not account for the periodic boundary
%     condition for sampling neurons.

if  is_simulation   
	Ic1=transpose(unique(s1(2,:)));
	Ic1=Ic1(Ic1<=Ne1^2);

else
	I1=transpose(unique(s1(2,:)));
	I1=I1(I1<=Ne1^2);
	Ix10=(ceil(I1/Ne1))/Ne1;
	Iy10=(mod((I1-1),Ne1)+1)/Ne1;
	I1=I1(Ix10<0.75 & Ix10>0.25 & Iy10<0.75 & Iy10>0.25);
    Ic1 = I1;
end