function ParamChange=configure_params(x,Ne1,Ni1,dt,T,is_small)
%% Generate the parameter struct for SNN simulation
%   -Input
%      x: table, parameter sets for the simulation
%      Ne1: int, number of E neurons in the recurrent layer per each dim
%      (total number of E neurons is thus Ne1^2)
%      Ni1: int, number of I neurons in the recurrent layer per each dim
%      dt: int, simulation time step length
%      T: int, simulation length
%      is_small: {0,1}, if using a smaller network size compared to (Huang
%      et al, 2019). Used for scaling the in-degree of network.
%   -Output
%      ParamChange: struct; configurations for the network simulation

rng shuffle;
Ntotal=Ne1*Ne1+Ni1*Ni1;
N_groundtruth=50000;
N_factor=sqrt(Ntotal/N_groundtruth); % May have to scale Jrs is using a smaller network


try 
	taudsynE=x.taudsynE;
catch 
	taudsynE=5;
end

try 
	taudsynI=x.taudsynI;
catch 
	taudsynI=8;
end

try 
	taursynE=x.taursynE;
catch 
	taursynE=1;
end

try 
	taursynI=x.taursynI;
catch 
	taursynI=1;
end

try 
	sigmaRREE=x.sigmaRREE;
catch 
	if is_small
		sigmaRREE=.1;
	else
		sigmaRREE=.1;
	end
end

try 
	sigmaRREI=x.sigmaRREI;
catch 
	if is_small
		sigmaRREI=.1;
	else
		sigmaRREI=.1;
	end
end

try 
	sigmaRRIE=x.sigmaRRIE;
catch 
	if is_small
		sigmaRRIE=.1;
	else
		sigmaRRIE=.1;
	end


end

try 
	sigmaRRII=x.sigmaRRII;
catch 
	
	if is_small
		sigmaRRII=.1;
	else
		sigmaRRII=.1;
	end

end

try 
	PrrEE=x.PrrEE;
catch 

	if is_small
		PrrEE=0.15;
	else
		PrrEE=.01;
	end

end

try 
	PrrEI=x.PrrEI;
catch 
	if is_small
		PrrEI=0.6;
	else
		PrrEI=.04;
	end

end

try 
	PrrIE=x.PrrIE;
catch 
	if is_small
		PrrIE= 0.45;
	else
		PrrIE=.03;
	end
end

try 
	PrrII=x.PrrII;
catch 

	if is_small
		PrrII=0.6;
	else
		PrrII=.04;
	end
end

try 
	JrEE=x.JrEE;
catch 
	JrEE=80*N_factor;
end
try 
	JrEI=x.JrEI;
catch 
	JrEI=-240*N_factor;
end
try 
	JrIE=x.JrIE;
catch 
	JrIE=40*N_factor;
end
try 
	JrII=x.JrII;
catch 
	JrII=-300*N_factor;
end
%%NEWLY ADDED
try 
	taudsynX=x.taudsynX;
catch 
	taudsynX=5;
end

try 
	taursynX=x.taursynX;
catch 
	taursynX=1;
end

try 
	sigmaREX=x.sigmaREX;
catch 
	sigmaREX=.05;
end

try 
	sigmaRIX=x.sigmaRIX;
catch 
	sigmaRIX=.05;
end

try 
	PsynX=x.PsynX;
catch 
	PsynX=1;
end

try 
	PsynE=x.PsynE;
catch 
	PsynE=1;
end

try 
	PsynI=x.PsynI;
catch 
	PsynI=1;
end

try 
	PrrEX=x.PrrEX;
catch 
	PrrEX=.1;
end

try 
	PrrIX=x.PrrIX;
catch 
	PrrIX=.05;
end

try 
	JrEX=x.JrEX;
catch
	if is_small
		JrEX=140*N_factor;
	else
		JrEX=140;
	end	 	
end

try 
	JrIX=x.JrIX;
catch 
	if is_small
		JrIX=100*N_factor;
	else
		JrIX=100;
	end	 	
end

try 
	PstimRX=x.PstimRX;
catch 
	PstimRX=.01;
end

try 
	Nx1=x.Nx1;
catch 
	Nx1=50;
end

try 
	Ne11=x.Ne11;
catch 
	Ne11=50;
end

try 
	Ni11=x.Ni11;
catch 
	Ni11=25;
end


try 
	mean_sigmaRRIs=x.mean_sigmaRRIs;
	sigmaRREI=mean_sigmaRRIs;
	sigmaRRII=mean_sigmaRRIs;

catch 
	1;
end


try 
	mean_sigmaRREs=x.mean_sigmaRREs;
	sigmaRREE=mean_sigmaRREs;
	sigmaRRIE=mean_sigmaRREs;

catch 
	1;
end

try 
	mean_sigmaRXs=x.mean_sigmaRXs;
	sigmaRIX=mean_sigmaRXs;
	sigmaREX=mean_sigmaRXs;

catch 
	1;
end


try 
	mean_JrIs=x.mean_JrIs;
	JrEI=mean_JrIs;
	JrII=mean_JrIs*1.25;

catch 
	1;
end


try 
	mean_PrrIs=x.mean_PrrIs;
	PrrEI=mean_PrrIs;
    PrrII=mean_PrrIs;
catch 
	1;
end

try 
	mean_JrXs=x.mean_JrXs;
	JrEX=mean_JrXs;
    JrIX=mean_JrXs*5/7;
catch 
	1;
end


filename='nothing';
ParamChange={'filename', filename;'dt', dt; 'T', T;'param(1).taudsyn',[taudsynX;taudsynE;taudsynI];'param(1).taursyn',[taursynX;taursynE;taursynI];'param(1).sigmaRR',[sigmaRREE,sigmaRREI;sigmaRRIE,sigmaRRII];'param(1).sigmaRX',[sigmaREX;sigmaRIX]; 'param(1).Prr',[PrrEE,PrrEI;PrrIE,PrrII]; 'param(1).Prx',[PrrEX;PrrIX]; 'param(1).Jr',[JrEE,JrEI;JrIE,JrII];'param(1).Jx',[JrEX;JrIX];'p_stim.rX',PstimRX;'param(1).Ne',Ne1*Ne1;'param(1).Ni',Ni1*Ni1;'param(1).N',Ntotal;'param(1).Nx',Nx1*Nx1};


