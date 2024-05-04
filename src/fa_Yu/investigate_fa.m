%apr3 2020

data_dim=10
latent_dim=2;
true_params.Ph = 0.01 * (1:data_dim)';
true_params.mu = ones(data_dim, 1);
true_params.d = ones(data_dim, 1);
L_random=randn(data_dim, latent_dim);
f = randn(data_dim,1);

true_params.L=[f,-f];

X = simdata_fa(true_params, 1000);

[estParams, ~]=fastfa(X,2);
[percentshared, d_shared, normevals] = compute_shared(estParams, 0.95);
COV=cov([X]');
Var=diag(COV);
R = COV./sqrt(Var*Var'); 
upper_R=R(triu(true(size(R)), 1));
mean_corr1=nanmean((upper_R));
d_shared
mean_corr1






data_dim=10
latent_dim=10;
true_params.Ph = 0.01 * (1:data_dim)';
true_params.mu = ones(data_dim, 1);
true_params.d = ones(data_dim, 1);
L_random=randn(data_dim, latent_dim);
f = randn(data_dim,1);

true_params.L=eye(10);

X = simdata_fa(true_params, 1000);

[estParams, ~]=fastfa(X,10);
[percentshared, d_shared, normevals] = compute_shared(estParams, 0.95);
COV=cov([X]');
Var=diag(COV);
R = COV./sqrt(Var*Var'); 
upper_R=R(triu(true(size(R)), 1));
mean_corr1=nanmean((upper_R));
d_shared
mean_corr1


figure
imshow(X(:,1:100))


figure
imshow(X_1(:,1:100))

X_low = X_1(:,1:100);

X_HIGH = X(:,1:100);

save('../results/sim_rsc_dim.mat','X_low','X_HIGH')
%just a little testing

data_dim=50
latent_dim=5;
true_params.Ph = 0.1 * (1:data_dim)';
true_params.mu = ones(data_dim, 1);
true_params.d = ones(data_dim, 1);
L_random=randn(data_dim, latent_dim);

Q=randn(data_dim, data_dim);
U=orth(Q);
%D=diag([39 14 1 1 1 zeros(1,data_dim-5)]);

D=diag([39 14 14 5 5 zeros(1,data_dim-latent_dim)]);

%D=diag([39 14 3 2 1 zeros(1,data_dim-5)]);

L=U*sqrt(D);

%L=L_sqrt*L_sqrt';

true_params.L=L;

X = simdata_fa(true_params, 350000);

for draw = 1:20
    X_train=X(:,((draw-1)*75+1):(draw*75));
    
    
    dim=crossvalidate_fa(X_train);
    istar = ([dim.sumLL] == max([dim.sumLL]));
    find((istar==1))
    %[estParams, LL]=fastfa(X_train_t,10);
    %[true_percentshared, true_d_shared,true_normevals] = compute_shared(estParams, 0.95);
    %true_d_shared

end





%%%%%%%%%%%%%%%%%%%%
%%%%%%%Feb%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%


data_dim=50；
latent_dim=5;
%true_params.L  = randn(50, 5);
true_params.Ph = 0.1 * (1:data_dim)';
true_params.mu = ones(data_dim, 1);
true_params.d = ones(data_dim, 1);


% create L for given eigen spectrum
Q=randn(data_dim, data_dim);
U=orth(Q);
D=diag([39 30 25 3 1 zeros(1,data_dim-latent_dim)]);
true_params.L=U*sqrt(D);
true_params.L=true_params.L(:,1:latent_dim);


X = simdata_fa(true_params, 600000);
X_validate = simdata_fa(true_params, 40000);
[real_eig_vectors, real_eigs]=eig(true_params.L*true_params.L');
[out,idx]=sort(diag(real_eigs),'descend');
real_eig_vectors=real_eig_vectors(:,idx);


[true_percentshared, true_d_shared,true_normevals] = compute_shared(true_params, 0.95);

%initialize statistics
sum_LLs = zeros(25,100,25);
eigen_spec=zeros(25,100);
d_shared_var=zeros(25,100);
percent_shared_var=zeros(25,100);
mu_diff=zeros(25,100);
psi_diff=zeros(25,100);
L_diff=zeros(25,100);
L_angles=zeros(25,100,5);
est_mus=zeros(25,100,data_dim);
est_psis=zeros(25,100,data_dim);
est_Ls=zeros(25,100,data_dim,5);



true_sum_LLs = zeros(25,1);


for draw=1:100
    draw
    X_train=X(:,((draw-1)*(2400+data_dim+1)+1):(draw*(2400+data_dim+1)));
    

    for t=1:25
        X_train_t = X_train(:,1:(data_dim+1+(t-1)*100));
        [estParams, LL]=fastfa(X_train_t,5);
        
        
        
        
        %cross-validating sum_LL
        for fold=1:5
            X_test=X_validate(:,((fold-1)*(2400+data_dim+1)+1):(fold*(2400+data_dim+1)));
            for tL=1:25
                X_test_tL=X_test(:,1:(data_dim+1+(tL-1)*100));
                [blah, LL] = fastfa_estep(X_test_tL, estParams);
                sum_LLs(t,draw,tL)=sum_LLs(t,draw,tL)+LL;
                if (draw==1) && (t==1) 
                    [blah, LL_true] = fastfa_estep(X_test_tL, true_params);

                    true_sum_LLs(tL)=true_sum_LLs(tL)+LL_true;
                end 
            end
             
        end
        
        for tL=1:25
            sum_LLs(t,draw,tL)=sum_LLs(t,draw,tL)/(5*(data_dim+1+(tL-1)*100));
            if (draw==1) && (t==1) 
             true_sum_LLs(tL)=true_sum_LLs(tL)/(5*(data_dim+1+(tL-1)*100));
            end 
        end
        
        %comparing fa dshared and eigenspectrum
        [percentshared, d_shared,normevals] = compute_shared(estParams, 0.95);
        normevalss_overlap=normevals(1:5);
        eigen_spec(t,draw)=sqrt(sum((normevals(1:5)/sum(normevals(1:5)) - true_normevals(1:5)/sum(true_normevals(1:5))).^2));
        eigen_spec_all(t,draw)=sqrt(sum((normevals(1:5)/sum(normevals(1:5)) - true_normevals(1:5)/sum(true_normevals(1:5))).^2));

        percent_shared_var(t,draw)=percentshared;  



        
        %comparing fa paras
        mu_diff(t,draw)=norm(estParams.d-true_params.d,'fro');
        psi_diff(t,draw)=norm(estParams.Ph-true_params.Ph,'fro');
        L_diff(t,draw)=norm(estParams.L-true_params.L,'fro');
        est_mus(t,draw,:)=estParams.d;
        est_psis(t,draw,:)=estParams.Ph;
        est_Ls(t,draw,:,:)=estParams.L;

        [est_eig_vectors, est_eigs]=eig(estParams.L*estParams.L');

        [out,idx]=sort(diag(est_eigs),'descend');
        est_eig_vectors=est_eig_vectors(:,idx);


        for ag = 1:5
        u = squeeze(est_eig_vectors(:,ag));    
        v = squeeze(real_eig_vectors(:,ag));
        L_angles(t,draw,ag)=min(abs(acosd(dot(u,v)/(norm(u)*norm(v)))) ,abs(acosd(dot(u,-v)/(norm(u)*norm(-v)))));
        end 

        %estimating dshared using m=10
        %[estParams, LL]=fastfa(X_train_t,10);
        %[percentshared, d_shared,normevals] = compute_shared(estParams, 0.95);

        %d_shared_var(t,draw)=d_shared; 

    end
end


fprintf('starting cv')

cv_m=zeros(80,100);
for draw = 1:100
    draw
    X_train=X(:,((draw-1)*5065+1):(draw*5065));
    for t=1:80
    X_train_t = X_train(:,1:(80+(t-1)*25));
    
    dim=crossvalidate_fa(X_train_t);
    istar = ([dim.sumLL] == max([dim.sumLL]));
    cv_m(t,draw) = find((istar==1));
    [estParams, LL]=fastfa(X_train_t,cv_m(t,draw));
    [percentshared, d_shared,normevals] = compute_shared(estParams, 0.95);
    d_shared_var(t,draw)=d_shared; 


    save('../investigate_fa_feb17_cvonly.mat','d_shared_var','true_d_shared', 'cv_m');


    end
end 



save('../investigate_fa_feb17_cvonly.mat','d_shared_var','true_d_shared', 'cv_m');

%screen 80217




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TEST latent vs real m%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dim=100；
latent_dim=10;
true_params.L  = randn(data_dim, latent_dim);
true_params.Ph = 0.1 * (1:data_dim)';
true_params.mu = ones(data_dim, 1);
true_params.d = ones(data_dim, 1);


% create L for given eigen spectrum
Q=randn(data_dim, data_dim);
U=orth(Q);
D=diag([30 25 20 15 10 1 1 1 1 1 zeros(1,data_dim-latent_dim)]);
true_params.L=U*sqrt(D);
true_params.L=true_params.L(:,1:latent_dim);



[real_eig_vectors, real_eigs]=eig(true_params.L*true_params.L');
[true_eig,idx]=sort(diag(real_eigs),'descend');
real_eig_vectors=real_eig_vectors(:,idx);


[true_percentshared, true_d_shared,true_normevals] = compute_shared(true_params, 0.95);




eigen_spec=zeros(16,100);
percent_of_shared=zeros(16,100);
percent_shared=zeros(16,100);
d_shared_var=zeros(16,100);
L_angles=zeros(16,100,5);



for draw=1:100
    draw
    Xtrain= simdata_fa(true_params, 20000);
    for m = 1:16
    
        [estParams, LL]=fastfa(Xtrain,m+4);


        [percentshared, d_shared,normevals] = compute_shared(estParams, 0.95);
       
        percent_shared(m,draw)=percentshared;  
        d_shared_var(m,draw)=d_shared;
        percent_of_shared(m,draw)=sqrt(sum((normevals(1:5) - true_normevals(1:5)).^2));
        

        


        [est_eig_vectors, est_eigs]=eig(estParams.L*estParams.L');

        [est_eigvals,idx]=sort(diag(est_eigs),'descend');
        est_eig_vectors=est_eig_vectors(:,idx);

        eigen_spec(m,draw)=sqrt(sum((est_eigvals(1:5) - true_eig(1:5)).^2));


        for ag = 1:5
        u = squeeze(est_eig_vectors(:,ag));    
        v = squeeze(real_eig_vectors(:,ag));
        L_angles(m,draw,ag)=min(abs(acosd(dot(u,v)/(norm(u)*norm(v)))) ,abs(acosd(dot(u,-v)/(norm(u)*norm(-v)))));
        end 

    end 

end

save('../investigate_fa_feb17_changem1.mat','eigen_spec','percent_of_shared','percent_shared','d_shared_var','L_angles','true_percentshared','true_d_shared');



latent_dim=10;


Q=randn(data_dim, data_dim);
U=orth(Q);
D=diag([550 500 450 400 350 250 250 250 250 250 zeros(1,data_dim-latent_dim)]);
true_params.L=U*sqrt(D);
true_params.L=true_params.L(:,1:latent_dim);



[real_eig_vectors, real_eigs]=eig(true_params.L*true_params.L');
[true_eig,idx]=sort(diag(real_eigs),'descend');
real_eig_vectors=real_eig_vectors(:,idx);


[true_percentshared, true_d_shared,true_normevals] = compute_shared(true_params, 0.95);




eigen_spec=zeros(16,100);
percent_of_shared=zeros(16,100);
percent_shared=zeros(16,100);
d_shared_var=zeros(16,100);
L_angles=zeros(16,100,5);


L_angles(:,1,1)
for draw=1:100
    draw
    Xtrain= simdata_fa(true_params, 20000);
    for m = 1:16
    
        [estParams, LL]=fastfa(Xtrain,m+4);


        [percentshared, d_shared,normevals] = compute_shared(estParams, 0.95);
       
        percent_shared(m,draw)=percentshared;  
        d_shared_var(m,draw)=d_shared;
        percent_of_shared(m,draw)=sqrt(sum((normevals(1:5) - true_normevals(1:5)).^2));
        

        


        [est_eig_vectors, est_eigs]=eig(estParams.L*estParams.L');

        [est_eigvals,idx]=sort(diag(est_eigs),'descend');
        est_eig_vectors=est_eig_vectors(:,idx);

        eigen_spec(m,draw)=sqrt(sum((est_eigvals(1:5) - true_eig(1:5)).^2));


        for ag = 1:5
        u = squeeze(est_eig_vectors(:,ag));    
        v = squeeze(real_eig_vectors(:,ag));
        L_angles(m,draw,ag)=min(abs(acosd(dot(u,v)/(norm(u)*norm(v)))) ,abs(acosd(dot(u,-v)/(norm(u)*norm(-v)))));
        end 

    end 

end

save('../investigate_fa_feb17_changem2.mat','eigen_spec','percent_of_shared','percent_shared','d_shared_var','L_angles','true_percentshared','true_d_shared');
%screen 26723


