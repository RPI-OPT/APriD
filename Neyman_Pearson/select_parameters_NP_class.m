%% code for Neyman-Pearson Classification Problem
% need a given filename before running
% need addpath ./data_set which include the data sets

%% 
switch filename
    case 'spambase'
        % algorithm parameter setting as stated in the paper
        opts.K = 1e5; f1_rp = 0.7; opts.maxit=5e5; opts.Jp = 10; opts.Jn = 10;
        opts.alpha_K = 10; opts.rho_K = 1; opts.gam_K= 10;
        ks_max_intervial=200;
%     case 'madelon'
%         % algorithm parameter setting as stated in the paper
%         opts.K = 1e5; f1_rp = 0.6; opts.maxit=5e5; opts.Jp = 10; opts.Jn = 10;
%         opts.alpha_K = 10; opts.rho_K = 1; opts.gam_K= 10;
%         ks_max_intervial=1000;
    case 'gisette'
        % this dataset has 5000 features, and and iALM  takes a long time. 
        % algorithm parameter setting as stated in the paper
        opts.K = 1e5; f1_rp = 0.6; opts.maxit=1e5; opts.Jp = 30; opts.Jn = 30;
        opts.alpha_K = 10; opts.rho_K = 1; opts.gam_K= 10;
        ks_max_intervial=2000;
    otherwise
        error('not the dataset or no the parameters for the dataset');
end
fprintf(['Neyman-Pearson Classification Problem on data set ' filename  '\n']);

% record the results at the iteration  k  in  ks
ks= 1; K=1e5;
for k=1+(1:K)
    if (mod(k,10^(floor(log10(k))))==0) ||  (mod(k,ks_max_intervial)==0) || k==K+1
        ks=[ks,k];
    end
end
opts.ks=ks;

opts.f1_r = -log(f1_rp); opts.eta = 0.04; 

%%
% load the data which is already normalized. 
load([filename '.mat']); 
M_p = X(y==1, :); M_n = X(y==-1,:);
m_p = size(M_p,1); m_n = size(M_n,1); d = size(M_p,2);
% initial with random x
opts.x = randn(1,d); 

%% iALM_NP
fprintf('iALM_NP is computing\n')
time = tic;
alpha = opts.f1_r; opts.w0 = (opts.x).'; 
[w, out] = iALM_NP(X,y,alpha,opts);
time_iALM = toc(time);
f0_opt = min(out.hist_obj(out.hist_res<=0));
f1_opt = min(out.hist_res);

if f1_opt>0
    fprintf('infeasible. \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Figure 4, theta = 10, alpha = 10, and rho in {100,10,1,0.1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% rho = 100
opts.theta = 10;  opts.alpha_K = 10; opts.rho_K = 100;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10_alpha10_rho100 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10_alpha10_rho100 = toc(time);
% rho = 10
opts.theta = 10;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10_alpha10_rho10 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10_alpha10_rho10 = toc(time);
% rho = 1 
opts.theta = 10;  opts.alpha_K = 10; opts.rho_K = 1;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10_alpha10_rho1 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10_alpha10_rho1 = toc(time); 
% rho = 0.1 
opts.theta = 10;  opts.alpha_K = 10; opts.rho_K = 0.1;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10_alpha10_rho10n1 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10_alpha10_rho10n1 = toc(time); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Figure 5, theta = 10, (alpha, rho) in {(10, 1), (5,2), (2,5), (1,10)}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% (alpha, rho) = (10,1) calculated in Figure 4 
% (alpha, rho) = (5,2)
opts.theta = 10; opts.alpha_K = 5; opts.rho_K = 2; 
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10_alpha5_rho2 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10_alpha5_rho2 = toc(time);
% (alpha, rho) = (2,5)
opts.theta = 10; opts.alpha_K = 2; opts.rho_K = 5; 
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10_alpha2_rho5 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10_alpha2_rho5 = toc(time);
% (alpha, rho) = (1,10)
opts.theta = 10; opts.alpha_K = 1; opts.rho_K = 10; 
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10_alpha1_rho10 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10_alpha1_rho10 = toc(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Figure 6, alpha=10 rho=10, theta in {1e-4,1e-3,1e-2,1e-1,1,1e1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% theta=1e-4
opts.theta = 0.0001;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10n4_alpha10_rho100 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10n4_alpha10_rho100 = toc(time);
% theta=1e-3
opts.theta = 0.001;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10n3_alpha10_rho100 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10n3_alpha10_rho100 = toc(time);
% theta=1e-2
opts.theta = 0.01;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10n2_alpha10_rho100 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10n2_alpha10_rho100 = toc(time);
% theta=1e-1
opts.theta = 0.1;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta10n1_alpha10_rho100 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta10n1_alpha10_rho100 = toc(time);
% theta=1
opts.theta = 1;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_NP_class_ApriD_theta1_alpha10_rho100 = NP_class_ApriD(M_p,M_n,opts);
time_NP_class_ApriD_theta1_alpha10_rho100 = toc(time);
% theta=10, calculated in Figure 4 

%%

%%
savename = [filename '_APriD_parameters'];
close all
clear X y M_n M_p time
save(savename) 

%% Plot Figures 4,5,6 in the paper with the above results.
plot_select_parameters_NP_class