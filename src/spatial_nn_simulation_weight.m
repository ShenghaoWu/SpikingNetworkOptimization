function [Wrr1,Wrf1,s1]=spatial_nn_simulation_weight(varargin)
%% Simulation script for the SBN/CBN. Authored by Chengcheng Huang

% RF2D3layer(option, ParamChange)

% param is a struc w/ fields: Ne, Ni, Nx, Jx, Jr, Kx, Kr, 
%       gl, Cm, vlb, vth, DeltaT, vT, vl, vre, tref, tausyn, V0, T, dt,
%       maxns, Irecord, Psyn
%   Jx=[Jex; Jix]; Jr=[Jee, Jei; Jie, Jii];
%   Kx=[Kex; Kix]; Kr=[Kee, Kei; Kie, Kii]; % out-degrees are fixed
%   taursyn: syn rise time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
%   taudsyn: syn decay time const, 3x(Nsyntype), rows: X, E, I; cols: syn type
%   Psyn(i,j): percentage of synapse j for (X, E, I (i=1,2,3))
%   sigmaRR=[sigmaee, sigmaei; sigmaie, sigmaii]; 
%   sigmaRX=[sigmaeX; sigmaiX]; 

%   Wrr is a vector of connections among the recurrent layer, containing postsynaptic cell indices, 
%       sorted by the index of the presynaptic cell. The block of postsynaptic cell indices for each presynaptic
%       cell is sorted as excitatory followed by inhibitory cells. I use fixed number of projections Kab to each population. 
%       For example, Wrr[j*(Kee+Kie)] to Wrr[j*{Kee+Kie)+Kee-1] are connections from j to E pop and 
%       Wrr[j*(Kee+Kie)+Kee] to Wrr[(j+1)*{Kee+Kie)-1] are connections from j to I pop. 
%   Wrf is a vector of connections from the feedforward layer to the recurrent layer, sorted by the index of the presynaptic cell.
%       The block of postsynaptic cell indices for each presynaptic cell is sorted as excitatory followed by inhibitory cells.

%   conversion of neuron ID (exc) to (x,y) coordinate in [1, Ne1]x[1, Ne1]: 
        % exc. ID [1, Ne], x=ceil(I/Ne1); y=(mod((I-1),Ne1)+1); ID=(x-1)*Ne1+y 
        % inh. ID [Ne+1, Ne+Ni], x=ceil((I-Ne)/Ni1); y=(mod((I-Ne-1),Ni1)+1); ID=(x-1)*Ni1+y+Ne; 
        
% sx: spike trains from Layer0
%     sx(1,:) contains spike times.
%     sx(2,:) contains indices of neurons that spike 
% s1: spike trains from Layer1
% s2: spike trains from Layer2
%  
% data save in filename 
% options is a struct w/ fields:
%   'save','CompCorr','plotPopR','fixW','loadS1','Layer1only'. Default values are 0.  
% ParamChange is a cell of 2 columns, 
%    the 1st column is variable names and
%    the 2nd column is the values. 

% if options.save is 1, ParamChange needs to have field 'filename'. 
% if options.CompCorr is 1, ParamChange needs to have field 'Nc',e.g. Nc=[500 500];
%      # of neurons to sample from Layer2 & Layer3. 
% if options.fixW is 1, ParamChange needs to have field 'Wseed1' & 'Wseed2'. 
rng shuffle;

nVarargs = length(varargin);
switch nVarargs
    case 1
        option = varargin{1};
    case 2
        option = varargin{1};
        ParamChange = varargin{2};
end 
    
if ~isfield(option, 'save') option.save=0; end 
if ~isfield(option, 'CompCorr') option.CompCorr=0; end 
if ~isfield(option, 'loadS1') option.loadS1=0; end 
if ~isfield(option, 'fixW') option.fixW=0; end 
if ~isfield(option, 'plotPopR') option.plotPopR=0; end 
if ~isfield(option, 'Layer1only') option.Layer1only=0; end 

if option.save==1 
    if ~ismember('filename',ParamChange(:,1))
        error('No filename to save data') 
    end 
end
if option.CompCorr==1 
    if ~ismember('Nc',ParamChange(:,1))
        error('No Nc (1x2): # of neurons to sample to compute correlations') 
    end  
end
if option.loadS1==1 
    if ~ismember('s1_fname',ParamChange(:,1))
        error('No s1_fname') 
    end 
end
if option.fixW==1 
    if ~ismember('Wseed1',ParamChange(:,1))|| ~ismember('Wseed2',ParamChange(:,1))
        error('No Wseed1 or Wseed2') 
    end   
end

%% define parameters 
dim ='2D';
% Number of neurons in network
Ne11=200;  % Number in each direction
Ni11=100; %100;
Ne21=200;  % Number in each direction
Ni21=100; %100;
Nx1=50;   % feedforward layer

param(1).Ne=Ne11*Ne11;
param(1).Ni=Ni11*Ni11;
param(1).Nx=Nx1*Nx1;
param(2).Ne=Ne21*Ne21;
param(2).Ni=Ni21*Ni21;
param(2).Nx=Ne11*Ne11;%feedforward from layer 2
param(1).N=param(1).Ne+param(1).Ni;
param(2).N=param(2).Ne+param(2).Ni;

% stimulus param
p_stim.Nstim=1;
p_stim.stim_type='Uncorr';
p_stim.rX=.01; % Rate of neurons in feedforward layer (kHz), size 1xNstim cell of element 1xNsource array
p_stim.Nsource=1;  % # of sources for global correlation, size 1xNstim
% p_stim.taucorr=10; % temporal jitter for 'LocalCorr' & 'GlobalCorr'
% p_stim.sigmac=0;  % local corr width
% p_stim.cx=0;    

% stim_type= 'LocalCorr'; correlation width 'sigmac'
% stim_type='spatialInput';  Gaussian inputs centered at 'center' with width 'sigmac' and mean rate 'rX' 
%   center=[.5 .5];
%   sigmac=.15;  % size 1xNstim
% stim_type= 'GlobalCorr'; % cell of size 1xNstim

T=20000; % Total sim time (in msec)

% static currents to Layer 3
inE=0;
inI=4;

% Connection widths
param(1).sigmaRX=.05*ones(2,1);
param(1).sigmaRR=.1*ones(2,2);%%%%%%%yes!
param(2).sigmaRX=.1*ones(2,1);
param(2).sigmaRR=.2*ones(2,2);

% number of neurons to record synaptic inputs and voltages from
nrecordE0=zeros(1,2); 
nrecordI0=zeros(1,2);

% Synaptic time constants 
param(2).taudsyn=[5 100; 5, 100; 8, 100]; % rows: X, E, I, column for different syn types
param(2).taursyn=[1 2; 1, 2; 1, 2]; % rows: X, E, I 
param(2).Psyn=[.2 .8; 1, 0; 1, 0]; % percentage of diff syn currents

param(1).taudsyn=[5; 5; 8]; % rows: X, E, I   %%%%%%%yes!
param(1).taursyn=[1; 1; 1]; %%%%%%%yes!
param(1).Psyn=[1; 1; 1];

% Connection probabilities (kind of, see use below)
param(1).Prr=[.01, .04; .03, .04];%%%%%%%yes!
param(2).Prr=[.01, .04; .03, .04];
param(1).Prx=[ .1; .05];
param(2).Prx=[ .05; .0];

% Connection strengths (scaled by sqrt(N) later)
%[Jee, Jei; Jie, Jii]
param(1).Jr=[80 -240; 40, -300];%%%%%%%yes!
param(2).Jr=[80 -240; 40, -300];
param(1).Jx=[140; 100]; 
param(2).Jx=[25; 0];

param(1).Iapp=[0;0];
param(2).Iapp=[inE;inI];

dt=.01;  % bin size  % 0.01
Tburn=1000;   % Burn-in period

% change parameters 
if nVarargs==2
    for i=1:size(ParamChange,1)
        eval([ParamChange{i,1} '= ParamChange{i,2};']);
    end
end
% param(2).Iapp=[inE;inI];
if option.verbose
    fprintf('\ninE=%.2f, inI=%.2f\n',param(2).Iapp(1),param(2).Iapp(2))
end
if option.loadS1
    p_stim.s1_fname=s1_fname;
end

%% initialization 

for par=1:2
    param(par).dt=dt;
    param(par).maxns=param(par).N*T*.06;
    param(par).T=T;
    % EIF neuron paramters
    param(par).gl=[1/15 1/10];  % E, I
    param(par).Cm=[1 1];
    param(par).vlb=[-100 -100];
    param(par).vth=[-10 -10];
    param(par).DeltaT=[2 .5];
    param(par).vT=[-50 -50]; %mV
    param(par).vre=[-65 -65];
    param(par).tref=[1.5 .5];
    V0min=param(par).vre(1);
    V0max=param(par).vT(1);
    param(par).vl=param(par).Iapp'.*[15, 10]-60;
    
    param(par).V0=(V0max-V0min).*rand(param(par).N,1)+V0min;
%     param(par).V0=data.param(par).V0;
    param(par).Kr=ceil(param(par).Prr.*[param(par).Ne, param(par).Ne; param(par).Ni,param(par).Ni]);
    param(par).Kx=ceil(param(par).Prx.*[param(par).Ne; param(par).Ni]);
    
    param(par).Irecord=[randi(param(par).Ne,1,nrecordE0(par)), (randi(param(par).Ni,1,nrecordI0(par))+param(par).Ne)];% neuron indice to record synaptic currents and Vm 
    param(par).Jr=param(par).Jr/sqrt(param(par).N);
    param(par).Jx=param(par).Jx/sqrt(param(par).N);
    
    % % Effective connection weights
    q=param(par).Ne/param(par).N;
    wrx=(param(par).Jx).*param(par).Prx*param(par).Nx/param(par).N;
    wrr=(param(par).Jr).*param(par).Prr.*[q, 1-q; q, 1-q];
    % For balanced state to exist this vector should be decreasing
    if option.verbose
        fprintf('\nThis list should be decreasing for\n  a balanced state to exist: %.2f %.2f %.2f\n\n',wrx(1)/wrx(2),abs(wrr(1,2)/wrr(2,2)),abs(wrr(1,1)/wrr(2,1)));
    end
    % and these values should be >1
    if option.verbose
        fprintf('\nAlso, this number should be greater than 1: %.2f\n\n',abs(wrr(2,2)/wrr(1,1)));
    end
    if par==1
        param(par).Imean=p_stim.rX*wrx*param(par).N+param(par).Iapp;
        if option.verbose
        disp(sprintf('\nFiring rates for large N: %.2f %.2f\n',-(wrr*param(par).N)\(p_stim.rX*wrx*param(par).N+param(par).Iapp)*1e3))
        end 
        firing_rates=-(wrr*param(par).N)\(p_stim.rX*wrx*param(par).N+param(par).Iapp)*1e3;
        

    end
    
end

%% generate input spike trains
if option.loadS1==0
   sx=genXspk(p_stim,param(1).Nx,T);
end

%% Simulation

% Simulate Network
if((param(2).N)<=200000)
    if option.verbose
        disp('simulation starts')
    end 
    % Random initial membrane potentials
    if option.loadS1
        load(s1_fname,'s1')
        s1=s1(:,s1(1,:)<=T);
        if option.verbose
            disp('load s1')
        end 
    else

        if option.givenW
            Ws=option.W1s;
            Wrr1=Ws{1};
            Wrf1=Ws{2};
        else

            if option.fixW
                param(1).Wseed=Wseed1;
                rng(Wseed1)
                if option.verbose
                    fprintf('seed%d for Wrr1, Wrf1\n',Wseed1)
                end
            end
            if option.verbose
                disp('generating Wrr1, Wrf1')
            end
            tic
            [Wrr1,Wrf1]=gen_weights(param(1).Ne,param(1).Ni,param(1).Nx,param(1).sigmaRX,param(1).sigmaRR,param(1).Prr, param(1).Prx,'2D'); 
            elapsetime=toc;
        end



        if option.verbose
            fprintf('elapsetime=%.2f sec\n',elapsetime)
        end
        if option.verbose
            disp('simulating Layer1')
        end
        tic
        [s1,Isyn1,Vm1]=EIF1DRFfastslowSyn(sx, Wrf1,Wrr1,param(1));
        s1=s1(:,s1(2,:)~=0);
        elapsetime=toc;
        if option.verbose
            fprintf('complete Layer1, elapsetime=%.2f sec\n',elapsetime)    
        end
    end
    if option.givenW
        Wrr1=0;
        Wrf1=0;
    end
else
    error('N too large') % Your computer probably can't handle this
    if option.verbose
        fprintf('N too large')
    end 
end

%save(filename,'s1','param','T')

end 


