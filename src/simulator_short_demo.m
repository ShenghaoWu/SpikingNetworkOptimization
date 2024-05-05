function [spike_train] =  simulator_short_demo(input_para,T)

%% Description: Function wrapper for the SNN simulator. This function wraps the actual simulation code, spatial_nn_simulation_weight.m, for the demonstration purpose in demo_short.m.
%% If you want to customize your own network model, wrap your simulator with a function of the exact same input/output formats of this function.
%% If you want to customize a different set of parameters for SBN, modify the parameter range and assign the named parameter with corresponding entry of input_para.

%   Input
%
%     input_para: [number of parameters to customize], parameter set to
%     simulate the network
%     T: int, network model simulation length (in ms).

%   Output
%     spike_train -- [2, number of spikes]. The first row is time (in ms), the second row is the neuron id. The first 2500 neurons are E neurons, the rest are I neurons.
%

taudsynI=input_para(1);
mean_sigmaRRIs = input_para(2);
% The rest of the parameters are set to the groundtruth values.
taudsynE = 22.102;
mean_sigmaRREs = 0.13862;
mean_sigmaRXs = 0.063017;
JrEI = -66.265;
JrIE = 27.965;
JrII = -74.641;
JrEE = 2.7119;
JrEX = 101.54;
JrIX = 35.767;
x= table(taudsynI,taudsynE,mean_sigmaRRIs,mean_sigmaRREs,mean_sigmaRXs,JrEI,JrIE,JrII,JrEE,JrEX,JrIX);

dt=0.05;
is_small=1;
Ne1=50;
Ni1=25;
ParamChange=configure_params(x,Ne1,Ni1,dt,T,is_small);

opt.save=0; % save data
opt.CompCorr=0; % compute correlations
opt.Layer1only=1; % 1 for two-layer network, 0 for three-layer network
opt.loadS1=0;
opt.plotPopR=0; % plot population rate
opt.fixW=0;  % use the same weight matrices for multiple simulations
opt.givenW=0;
opt.verbose=0;
opt.check_firing_rates=0;


[~,~,spike_train]=spatial_nn_simulation_weight(opt, ParamChange);


