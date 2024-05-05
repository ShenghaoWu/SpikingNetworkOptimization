function SNOPS_short_demo(target_stats_name,simulator,parameter_range,max_iter,T,save_name,is_plot,varargin)
%% SNOPS short demo function. This script is similar to bayesian_optimization.m but incorporates printing and plotting functionality for the demonstration purpose in demo_short.m.
%
%   -Input
%     target_stats_name: string; filename of the target activity statistics to customize network model to.  
%     simulator: function; the network model simulator.
%     parameter_range: [number of parameters, 2], Lower and upper
%           bounds of the search range for the parameters.
%     max_iter: int, number of iterations to run SNOPS.
%     T: int, network model simulation length (in ms).
%     save_name: string, filename for the log files.
%     is_plot: {0,1}, if =1, plot the customization info. Otherwise
%           the info will be printed to the console.


warning('off')
addpath('./src')
addpath('./src/fa_Yu')


[obj_configs,optimization_opt]=generate_config('T', T,...
                                               'max_iter',max_iter,...
                                               'x_range',parameter_range,...
                                               'min_cost_eval',1,...
                                               'max_cost_eval',1,...
                                               'real_data_name', target_stats_name,...
                                               'filename',save_name);
                                         
                                          

rng(42);
tt=obj_configs.true_statistics;
disp('Customization begins')
disp('================================================================================')
disp('Target activity statistics:')
disp(strcat('firing rate: ', string(tt.rate_mean),' || Fano: ', string(tt.fano_mean),' || rsc: ', string(tt.mean_corr_mean)))
disp(strcat('percent shared: ', string(100*tt.fa_percent_mean),' || d shared: ', string(tt.fa_dim_mean),' || eigenspectrum(top3): ', strjoin(string(tt.fa_normeval_mean(1:3)),',')))
pause(1)

%preparing cost function
func=@(x,is_surrogate)cost_func_short_demo(x,is_surrogate,obj_configs,simulator);

%initialization
x_range=optimization_opt.x_range;
epsilon=optimization_opt.epsilon;
n_sample=optimization_opt.n_sample;
n_local=optimization_opt.n_local;
max_time=optimization_opt.max_time;
max_iter=optimization_opt.max_iter-optimization_opt.n_data;
save_name=optimization_opt.save_name;
std_tol=optimization_opt.std_tol;
max_cost_eval=optimization_opt.max_cost_eval;
min_cost_eval=optimization_opt.min_cost_eval;
is_log= optimization_opt.is_log;
n_dim=size(x_range,1);


try 
	x_train=optimization_opt.x_train;
	y_train=optimization_opt.y_train;
	incumbent_std=optimization_opt.incumbent_std;
	n_data=size(x_train,1);
	optimization_time=optimization_opt.optimization_time;
	y_feasibility=optimization_opt.y_feasibility;

catch
   
	incumbent_std=0;
	n_data=optimization_opt.n_data;
	optimization_time=zeros(n_data,1);
	y_feasibility=NaN(n_data,1);
	try 
	    x_init=optimization_opt.x_init;  
    catch
        x_init=rand(n_data,n_dim).*repmat((x_range(:,2)-x_range(:,1))', n_data,1)+ repmat(x_range(:,1)', n_data,1);
    end
    
	%evaluating initial parameter sets
	x_train=x_init;
	y_train=NaN(n_data,1);
    y_train_best = NaN(n_data,1);
    
	for i =1:n_data
	    t_init=tic;
	    y_hat=func(x_train(i,:),1);
	    if isnan(y_hat)
	        y_feasibility(i)=0;
	    else 
	    	if isnan(nanmin(y_train)) || nanmin(y_train)+2*incumbent_std>y_hat
		    	intens_cnt = 1;
		    	tmp_yhat=[y_hat];
		    	current_std=Inf;
                while (intens_cnt< min_cost_eval) && (~any(isnan(tmp_yhat)))&&((mean(tmp_yhat)<nanmin(y_train)+2*incumbent_std)||isnan(nanmin(y_train)))
		    		intens_cnt=intens_cnt+1;
		    		tmp_yhat=[tmp_yhat,func(x_train(i,:),1)];
		    		current_std=std(tmp_yhat)/sqrt(intens_cnt);
		    	end
		    	while (intens_cnt< max_cost_eval) && current_std>std_tol && (~any(isnan(tmp_yhat)))&&((mean(tmp_yhat)<nanmin(y_train)+2*incumbent_std)||isnan(nanmin(y_train)))
		    		intens_cnt=intens_cnt+1;
		    		tmp_yhat=[tmp_yhat,func(x_train(i,:),1)];
		    		current_std=std(tmp_yhat)/sqrt(intens_cnt);
		    	end
		    	y_hat=mean(tmp_yhat);
		    	if y_hat < nanmin(y_train) || (isnan(nanmin(y_train))&& ~isnan(y_hat))
		    		incumbent_std=current_std;
		    	end
			end
			if isnan(y_hat)
		        y_feasibility(i)=0;
		    else
		        y_feasibility(i)=1;
		        y_train(i)=y_hat;
		    end
	    end
	    optimization_time(i)=toc(t_init);

    
    
    
    y_train_best(i)=nanmin(y_train);
    ts = obj_configs.true_statistics;
    load(strcat('./results/',obj_configs.filename,'_stats.mat'))  
    [I,J] = min(y_train);
    tt = full_stats(J,:);
    
    if is_plot
        try
            close(ff)
        catch
        end
        ff = figure('Name',strcat('iteration', {' '},string(i)),'NumberTitle','off');
        set(gcf, 'Position',  [500, 400, 1200, 1000])
        tiledlayout(2,8)
        
        ax1 = nexttile;
        hold(ax1,'on')
        scatter(0.01,tt.rate1,35,'black','filled')
        scatter(0,ts.rate_mean,35,'red','filled')
        xlabel('$fr$','Interpreter','latex')
        set(ax1,'XTick',[])
        hold(ax1,'off')

        ax2 = nexttile;
        hold(ax2,'on')
        scatter(0.01,tt.FanoFactor1,35,'black','filled')
        scatter(0,ts.fano_mean,35,'red','filled')
        xlabel('$ff$','Interpreter','latex')
        set(ax2,'XTick',[])
        hold(ax2,'off')

        ax3 = nexttile;
        hold(ax3,'on')
        scatter(0.01,tt.mean_corr1,35,'black','filled')
        scatter(0,ts.mean_corr_mean,35,'red','filled')
        xlabel('$r_{sc}$','Interpreter','latex')
        set(ax3,'XTick',[])
        hold(ax3,'off')

        ax4 = nexttile;
        hold(ax4,'on')
        scatter(0.01,100*tt.fa_percentshared,35,'black','filled')
        scatter(0,100*ts.fa_percent_mean,35,'red','filled')
        xlabel('$\%_{sh}$','Interpreter','latex')
        set(ax4,'XTick',[])
        hold(ax4,'off')

        ax5 = nexttile;
        hold(ax5,'on')
        scatter(0.01,tt.fa_dshared,35,'black','filled')
        scatter(0,ts.fa_dim_mean,35,'red','filled')
        xlabel('$d_{sh}$','Interpreter','latex')
        set(ax5,'XTick',[])
        hold(ax5,'off')

        ax6 = nexttile([1,3]);
        hold(ax6,'on')
        scatter(0.01+(1:3),tt.fa_normevals(1:3),35,'black','filled')
        scatter(1:3,ts.fa_normeval_mean(1:3),35,'red','filled')
        plot(0.01+(1:3),tt.fa_normevals(1:3),'LineWidth',2.5,'Color','black')
        plot(1:3,ts.fa_normeval_mean(1:3),'LineWidth',2.5,'Color','red')
        xlabel('$es$ (Top 3 modes)' ,'Interpreter','latex')
        legend('best','target','Location','northeast')
        set(ax6,'XTick',[])
        hold(ax6,'off')
        sgtitle('Activity statistics','FontSize',20,'FontWeight','bold'); 

        
        ax7 = nexttile([1,3]);
        hold(ax7,'on')
        scatter(x_train(1:i,1),x_train(1:i,2),35,[0 0 0]+0.05*12,'filled')
        scatter(x_train(i,1),x_train(i,2),35,'black','filled')
        title('Evaluated parameter sets')
        xlabel('$\tau^{id}$','Interpreter','latex')
        ylabel('$\sigma^{i}$','Interpreter','latex')
        xlim([0,35])
        ylim([0,0.35])
        hold(ax7,'off')

        
        ax8 = nexttile([1,5]);
        plot(1:length(y_train_best),log(y_train_best),'LineWidth',2.5,'Color','black');
        hold(ax8,'on')
        scatter(1:length(y_train_best),log(y_train_best),35,'black','filled'); 
        xlim([1,optimization_opt.max_iter+1])
        title('Cost trace')
        xlabel('iteration')
        ylabel('log best(cost)')
        hold(ax8,'off') 
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        drawnow

    end
    

    disp('================================================================================')
    if y_feasibility(i)==1
        load(strcat('./results/',obj_configs.filename,'_stats.mat'))  
        tt = full_stats(end,:);
        disp(strcat('iteration ', {' '},string(i)))
        disp('--------------------------------------------------------------------------------')
        disp(strcat('current parameter set:  ',join(string(x_train(i,:)),',')))
        disp('--------------------------------------------------------------------------------')
        disp(strcat('current cost (log)=',string(log(y_train(i))),'; || current best (log)=',string(log(nanmin(y_train)))))
        disp('--------------------------------------------------------------------------------')
        disp(strcat('firing rate: ', string(tt.rate1),' || Fano: ', string(tt.FanoFactor1),' || rsc: ', string(tt.mean_corr1)))
        disp(strcat('percent shared: ', string(100*tt.fa_percentshared),' || dshared: ', string(tt.fa_dshared),' || eigenspectrum(top3): ', strjoin(string(tt.fa_normevals(1:3)),',')))

    else 
        disp(strcat('iteration ', {' '}, string(i)))
        disp('--------------------------------------------------------------------------------')
        disp(strcat('current parameter set:  ',join(string(x_train(i,:)),',')))
        disp('--------------------------------------------------------------------------------')
        disp( 'infeasible parameter set')
    end
       
end

end 






%in case none of the initial points are feasible, gprMdl posterior is uniform 
if sum(y_feasibility==0)==n_data
    current_feas=0;
    gprMdl = fitrgp(x_train(1,:),0,...
                'KernelFunction','ardmatern52',...
                'BasisFunction', 'constant', ...
                'Standardize', false,...
                'ConstantSigma',false,...
                'FitMethod','exact','PredictMethod','exact');
        f_plus=0;
        epsilon_scaled=epsilon;
else
    current_feas=1;
	gprMdl=nan;
	f_plus=nan;
	epsilon_scaled=nan;

end

save(save_name,'x_train','y_train','y_feasibility','optimization_time','incumbent_std')


%main BO loop
for kk=1:max_iter-1
    if sum(optimization_time)<=max_time
        t_init=tic;
        
        [gprMdl,gprMdl_feas,f_plus,epsilon_scaled]=fit_gp(x_train,y_train,y_feasibility,epsilon,gprMdl,current_feas,f_plus,epsilon_scaled,is_log);
        [EI_maxx,~] = max_EI(gprMdl,gprMdl_feas,f_plus,epsilon_scaled,n_sample,x_range,n_local);
        x_train(n_data+kk,:)=EI_maxx;
        y_hat=func(EI_maxx,1);
        
        if isnan(y_hat)
            y_feasibility(n_data+kk)=0;
            y_train(n_data+kk)=y_hat;
            current_feas=0;
        else
	        if isnan(nanmin(y_train)) || nanmin(y_train)+2*incumbent_std>y_hat
		    	intens_cnt = 1;
		    	tmp_yhat=[y_hat];
		    	current_std=Inf;
                while (intens_cnt< min_cost_eval) && (~any(isnan(tmp_yhat)))&&((mean(tmp_yhat)<nanmin(y_train)+2*incumbent_std)||isnan(nanmin(y_train))) 
		    		intens_cnt=intens_cnt+1;
		    		tmp_yhat=[tmp_yhat,func(x_train(n_data+kk,:),1)];
		    		current_std=std(tmp_yhat)/sqrt(intens_cnt);
		    	end
		    	while (intens_cnt< max_cost_eval) && current_std>std_tol && (~any(isnan(tmp_yhat)))&&((mean(tmp_yhat)<nanmin(y_train)+2*incumbent_std)||isnan(nanmin(y_train))) 
		    		intens_cnt=intens_cnt+1;
		    		tmp_yhat=[tmp_yhat,func(x_train(n_data+kk,:),1)];
		    		current_std=std(tmp_yhat)/sqrt(intens_cnt);
		    	end
		    	y_hat=mean(tmp_yhat);
		    	if y_hat < nanmin(y_train)|| (isnan(nanmin(y_train))&& ~isnan(y_hat))
		    		incumbent_std=current_std;
		    	end
			end
			if isnan(y_hat)
		        y_feasibility(n_data+kk)=0;
		        y_train(n_data+kk)=y_hat;
		        current_feas=0;
		    else
		        y_feasibility(n_data+kk)=1;
		        y_train(n_data+kk)=y_hat;
		        current_feas=1;
		    end

        end
        optimization_time(n_data+kk)=toc(t_init);
        save(save_name,'x_train','y_train','y_feasibility','optimization_time','incumbent_std')
    else
        save(save_name,'x_train','y_train','y_feasibility','optimization_time','incumbent_std')
        return
    end
   

    y_train_best(end+1)=nanmin(y_train);
    load(strcat('./data/',obj_configs.filename,'_stats.mat'))  
    [I,J] = min(y_train);
    tt = full_stats(J,:);
   
    if is_plot
        try
            close(ff)
        catch
        end
        ff = figure('Name',strcat('iteration ', {' '},string(size(x_train,1))),'NumberTitle','off');
        set(gcf, 'Position',  [500, 400, 1200, 1000])
        tiledlayout(2,8)

        ax1 = nexttile;
        hold(ax1,'on')
        scatter(0.01,tt.rate1,35,'black','filled')
        scatter(0,ts.rate_mean,35,'red','filled')
        xlabel('$fr$','Interpreter','latex')
        set(ax1,'XTick',[])
        hold(ax1,'off')

        ax2 = nexttile;
        hold(ax2,'on')
        scatter(0.01,tt.FanoFactor1,35,'black','filled')
        scatter(0,ts.fano_mean,35,'red','filled')
        xlabel('$ff$','Interpreter','latex')
        set(ax2,'XTick',[])
        hold(ax2,'off')

        ax3 = nexttile;
        hold(ax3,'on')
        scatter(0.01,tt.mean_corr1,35,'black','filled')
        scatter(0,ts.mean_corr_mean,35,'red','filled')
        xlabel('$r_{sc}$','Interpreter','latex')
        set(ax3,'XTick',[])
        hold(ax3,'off')

        ax4 = nexttile;
        hold(ax4,'on')
        scatter(0.01,100*tt.fa_percentshared,35,'black','filled')
        scatter(0,100*ts.fa_percent_mean,35,'red','filled')
        xlabel('$\%_{sh}$','Interpreter','latex')
        set(ax4,'XTick',[])
        hold(ax4,'off')

        ax5 = nexttile;
        hold(ax5,'on')
        scatter(0.01,tt.fa_dshared,35,'black','filled')
        scatter(0,ts.fa_dim_mean,35,'red','filled')
        xlabel('$d_{sh}$','Interpreter','latex')
        set(ax5,'XTick',[])
        hold(ax5,'off')

        ax6 = nexttile([1,3]);
        hold(ax6,'on')
        scatter(0.01+(1:3),tt.fa_normevals(1:3),35,'black','filled')
        scatter(1:3,ts.fa_normeval_mean(1:3),35,'red','filled')
        plot(0.01+(1:3),tt.fa_normevals(1:3),'LineWidth',2.5,'Color','black')
        plot(1:3,ts.fa_normeval_mean(1:3),'LineWidth',2.5,'Color','red')
        xlabel('$es$ (Top 3 modes)' ,'Interpreter','latex')
        legend('best','target','Location','northeast')
        set(ax6,'XTick',[])
        hold(ax6,'off')
        sgtitle('Activity statistics','FontSize',20,'FontWeight','bold'); 
        
        ax7 = nexttile([1,3]);
        hold(ax7,'on')
        scatter(x_train(:,1),x_train(:,2),35,[0 0 0]+0.05*12,'filled')
        scatter(x_train(end,1),x_train(end,2),35,'black','filled')
        title('Evaluated parameter sets')
        xlabel('$\tau^{id}$','Interpreter','latex')
        ylabel('$\sigma^{i}$','Interpreter','latex')
        xlim([0,35])
        ylim([0,0.35])
        hold(ax7,'off')

        ax8 = nexttile([1,5]);
        plot(1:length(y_train_best),log(y_train_best),'LineWidth',2.5,'Color','black');
        hold(ax8,'on')
        scatter(1:length(y_train_best),log(y_train_best),35,'black','filled');  
        xlim([1,optimization_opt.max_iter+1])
        title('Cost trace')
        xlabel('iteration')
        ylabel('log best(cost)')
        hold(ax8,'off') 
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        drawnow
       
    end 
    disp('================================================================================')
    if y_feasibility(end)==1
        load(strcat('./data/',obj_configs.filename,'_stats.mat'))  
        tt = full_stats(end,:);
        disp(strcat('iteration ', {' '}, string(size(x_train,1))))
        disp('--------------------------------------------------------------------------------')
        disp(strcat('current parameter set:  ',join(string(x_train(end,:)),',')))
        disp('--------------------------------------------------------------------------------')
        disp(strcat('current cost (log)=',string(log(y_train(end))),'; || current best (log)=',string(log(nanmin(y_train)))))
        disp('--------------------------------------------------------------------------------')
        disp(strcat('firing rate: ', string(tt.rate1),' || Fano: ', string(tt.FanoFactor1),' || rsc: ', string(tt.mean_corr1)))
        disp(strcat('percent shared: ', string(100*tt.fa_percentshared),' || d shared: ', string(tt.fa_dshared),' || eigenspectrum(top3): ', strjoin(string(tt.fa_normevals(1:3)),',')))

    else 
        disp(strcat('iteration ', {' '}, string(size(x_train,1))))
        disp('--------------------------------------------------------------------------------')
        disp(strcat('current parameter set:  ',join(string(x_train(end,:)),',')))
        disp('--------------------------------------------------------------------------------')
        disp( 'infeasible parameter set')
    end
    
end

load(strcat('./data/',obj_configs.filename,'_stats.mat'))  
[I,J] = min(y_train);
tt = full_stats(J,:);
ps = paras(J,:);

 
disp('================================================================================')
disp(strcat('Customization complete. Total time elapsed (s): ', string(sum(optimization_time))))
disp('--------------------------------------------------------------------------------')
disp(strcat('best parameter set:  ',join(string(x_train(J,:)),',')))
disp('--------------------------------------------------------------------------------')
disp(strcat('best cost =',string(log(nanmin(y_train)))))
disp('--------------------------------------------------------------------------------')
disp(strcat('firing rate: ', string(tt.rate1),' || Fano: ', string(tt.FanoFactor1),' || rsc: ', string(tt.mean_corr1)))
disp(strcat('percent shared: ', string(100*tt.fa_percentshared),' || d shared: ', string(tt.fa_dshared),' || eigenspectrum(top3): ', strjoin(string(tt.fa_normevals(1:3)),',')))

if is_plot
        pause(3)
        try
            close(ff)
        catch
        end
        ff = figure('Name','Customization complete','NumberTitle','off');
        set(gcf, 'Position',  [500, 400, 1200, 1000])
        tiledlayout(2,8)

        ax1 = nexttile;
        hold(ax1,'on')
        scatter(0.01,tt.rate1,35,'black','filled')
        scatter(0,ts.rate_mean,35,'red','filled')
        xlabel('$fr$','Interpreter','latex')
        ylim([0,10])
        set(ax1,'XTick',[])
        hold(ax1,'off')

        ax2 = nexttile;
        hold(ax2,'on')
        scatter(0.01,tt.FanoFactor1,35,'black','filled')
        scatter(0,ts.fano_mean,35,'red','filled')
        xlabel('$ff$','Interpreter','latex')
        ylim([0,1.2])
        set(ax2,'XTick',[])
        hold(ax2,'off')

        ax3 = nexttile;
        hold(ax3,'on')
        scatter(0.01,tt.mean_corr1,35,'black','filled')
        scatter(0,ts.mean_corr_mean,35,'red','filled')
        xlabel('$r_{sc}$','Interpreter','latex')
        ylim([-0.05,0.05])
        set(ax3,'XTick',[])
        hold(ax3,'off')

        ax4 = nexttile;
        hold(ax4,'on')
        scatter(0.01,100*tt.fa_percentshared,35,'black','filled')
        scatter(0,100*ts.fa_percent_mean,35,'red','filled')
        xlabel('$\%_{sh}$','Interpreter','latex')
        ylim([0,100])
        set(ax4,'XTick',[])
        hold(ax4,'off')

        ax5 = nexttile;
        hold(ax5,'on')
        scatter(0.01,tt.fa_dshared,35,'black','filled')
        scatter(0,ts.fa_dim_mean,35,'red','filled')
        xlabel('$d_{sh}$','Interpreter','latex')
        set(ax5,'XTick',[])
        hold(ax5,'off')

        ax6 = nexttile([1,3]);
        hold(ax6,'on')
        scatter(0.01+(1:3),tt.fa_normevals(1:3),35,'black','filled')
        scatter(1:3,ts.fa_normeval_mean(1:3),35,'red','filled')
        plot(0.01+(1:3),tt.fa_normevals(1:3),'LineWidth',2.5,'Color','black')
        plot(1:3,ts.fa_normeval_mean(1:3),'LineWidth',2.5,'Color','red')
        xlabel('$es$ (Top 3 modes)' ,'Interpreter','latex')
        legend('best','target','Location','northeast')
        set(ax6,'XTick',[])
        hold(ax6,'off')
        sgtitle('Activity statistics','FontSize',20,'FontWeight','bold'); 
        
        ax7 = nexttile([1,3]);
        hold(ax7,'on')
        scatter(x_train(:,1),x_train(:,2),35,[0 0 0]+0.05*12,'filled')
        scatter(x_train(J,1),x_train(J,2),35,'black','filled')
        title('Evaluated parameter sets')
        xlabel('$\tau^{id}$','Interpreter','latex')
        ylabel('$\sigma^{i}$','Interpreter','latex')
        xlim([0,35])
        ylim([0,0.35])
        hold(ax7,'off')

        ax8 = nexttile([1,5]);
        plot(1:length(y_train_best),log(y_train_best),'LineWidth',2.5,'Color','black');
        hold(ax8,'on')
        scatter(1:length(y_train_best),log(y_train_best),35,'black','filled');  
        xlim([1,optimization_opt.max_iter+1])
        title('Cost trace')
        xlabel('iteration')
        ylabel('log best(cost)')
        hold(ax8,'off') 
        set(findall(gcf,'-property','FontSize'),'FontSize',20)
        drawnow
       
end 
    

end