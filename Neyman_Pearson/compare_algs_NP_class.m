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
    case 'madelon'
        % algorithm parameter setting as stated in the paper
        opts.K = 1e5; f1_rp = 0.6; opts.maxit=5e5; opts.Jp = 10; opts.Jn = 10;
        opts.alpha_K = 10; opts.rho_K = 1; opts.gam_K= 10;
        ks_max_intervial=1000;
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

% record the results at the iteration k in ks
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
%%
opts.x = randn(1,d); 

%% calculate the exact solution by iALM algorithm
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
%% run the ApriD algorithm
fprintf('ApriD is computing\n'); 
time = tic;
out_ApriD = NP_class_ApriD(M_p,M_n,opts);
time_ApriD = toc(time);
%%%%%
%% run the CSA algorithm
fprintf('CSA is computing\n'); 
time = tic;
out_CSA = NP_class_CSA(M_p,M_n,opts);
time_CSA = toc(time); 
%%%%%
%% run the MSA algorithm
% opts.rho_K  = opts.alpha_K;
fprintf('MSA is computing\n'); 
time = tic;
out_MSA = NP_class_MSA(M_p,M_n,opts);
time_MSA = toc(time); 
%%%%%
%% save the results 
savename = [filename '_ks_K_' num2str(opts.K) '_alpha_' num2str(floor(opts.alpha_K)) ...
    '_10rho_' num2str(floor(10*opts.rho_K)) '_gam_' num2str(opts.gam_K)  '_f1right_10prob'  num2str(10*f1_rp)  '_1000eta_' num2str(1000*opts.eta)];
%
close all
clear X y M_n M_p time
save(savename)
%% print and save the running time
% print the running time
fprintf([filename ' iALM takes time ' num2str(time_iALM) '\n']) 
 
fprintf(['computing time with ' filename  '\n'])
fprintf(['Example' ' & ' 'data or size' ' & ' 'APriD' ' & ' 'MSA' ' & ' 'CSA' ' & '  'PDSG_adp'  '  \\\\\n'])
fprintf(['NPC(5.2)' ' & ' filename ' & ' num2str(time_ApriD,'%.1f') ' & ' num2str(time_CSA,'%.1f') ' & ' num2str(time_MSA,'%.1f') ' & ' '-' '  \\\\\n'])

% save the running time of this example to the file running_time.txt.
file_time = fopen('running_time.txt','a');
fprintf(file_time, ['NPC(5.2)' ' & ' filename ' & ' num2str(time_ApriD,'%.1f') ' & ' num2str(time_CSA,'%.1f') ' & ' num2str(time_MSA,'%.1f') ' & ' '-' '  \\\\\n']);
fclose(file_time);

%% Plot Figures 1 in the paper with the above results.
plot_compare_algs_NP_class 