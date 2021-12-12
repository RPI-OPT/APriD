%% code for QCQP in Expectation Form
%  need a given n before  running

%%% min_{x in [X_min,X_max]^n} E[1/2\|H_i*x-c_i\|^2]
%%% s.t.  E[1/2 x^TQ_jx+a_j^Tx-b_j]<=0, 

%%
switch n
    case 10 
        p=5; % n = 10; p = 5;
        opts.K = 50000;  ks_max_intervial = 10; 
        fprintf(['QCQP in Expectation Form with n=' num2str(n) ', p=' num2str(p) '\n']);
    case 200
        p=150; %n = 200; p = 150;
        opts.K = 50000;  ks_max_intervial = 10;
        fprintf(['QCQP in Expectation Form with n=' num2str(n) ', p=' num2str(p) '\n']);
    otherwise
        error('new n, need give a new p');
end

N = 1e5;  M = 1e5;
X_min = -10; X_max = 10;

% call for data generate 
[EHTH,EcTH,EcTc,EQ,Ea,Eb] = data_generate_QCQP_exp(n,p,N,M); 

%%
opts.x = randn(1,n);
while 0.5*(opts.x)*EQ*(opts.x).'+(opts.x)*Ea-Eb<=1.0 %whether feasible
    opts.x = randn(1,n);
end
 
% record the results at the iteration  k  in  ks
% in  this  example, recording does take much time.
ks= 1;
for k=1+(1:opts.K)
    if (mod(k,10^(floor(log10(k))))==0) ||  (mod(k,ks_max_intervial)==0) || k==opts.K+1
        ks=[ks,k];
    end
end
opts.ks=ks;

%%
fprintf('data generated; begin computing. \n')
fprintf('cvx  is computing \n')
time = tic;
cvx_solver Mosek
cvx_begin
variable x(n)
minimize(0.5*(quad_form(x,EHTH)-2*EcTH*x+EcTc));
0.5*quad_form(x,EQ)+Ea.'*x-Eb <= 0;
X_min <= x <= X_max
cvx_end
time_expect_cvx = toc(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Figure 4, theta = 10, alpha = 10, and rho in {100,10,1,0.1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% rho = 100
opts.theta = 10;  opts.alpha_K = 10; opts.rho_K = 100;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10_alpha10_rho100 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10_alpha10_rho100 = toc(time);
% rho = 10
opts.theta = 10;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10_alpha10_rho10 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10_alpha10_rho10 = toc(time);
% rho = 1 
opts.theta = 10;  opts.alpha_K = 10; opts.rho_K = 1;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10_alpha10_rho1 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10_alpha10_rho1 = toc(time); 
% rho = 0.1 
opts.theta = 10;  opts.alpha_K = 10; opts.rho_K = 0.1;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10_alpha10_rho10n1 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10_alpha10_rho10n1 = toc(time); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Figure 5, theta = 10, (alpha, rho) in {(10, 1), (5,2), (2,5), (1,10)}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% (alpha, rho) = (10,1) calculated in Figure 4 
% (alpha, rho) = (5,2)
opts.theta = 10; opts.alpha_K = 5; opts.rho_K = 2; 
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10_alpha5_rho2 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10_alpha5_rho2 = toc(time);
% (alpha, rho) = (2,5)
opts.theta = 10; opts.alpha_K = 2; opts.rho_K = 5; 
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10_alpha2_rho5 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10_alpha2_rho5 = toc(time);
% (alpha, rho) = (1,10)
opts.theta = 10; opts.alpha_K = 1; opts.rho_K = 10; 
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10_alpha1_rho10 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10_alpha1_rho10 = toc(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Figure 6, alpha=10 rho=10, theta in {1e-4,1e-3,1e-2,1e-1,1,1e1}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% theta=1e-4
opts.theta = 0.0001;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10n4_alpha10_rho100 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10n4_alpha10_rho100 = toc(time);
% theta=1e-3
opts.theta = 0.001;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10n3_alpha10_rho100 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10n3_alpha10_rho100 = toc(time);
% theta=1e-2
opts.theta = 0.01;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10n2_alpha10_rho100 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10n2_alpha10_rho100 = toc(time);
% theta=1e-1
opts.theta = 0.1;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta10n1_alpha10_rho100 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta10n1_alpha10_rho100 = toc(time);
% theta=1
opts.theta = 1;  opts.alpha_K = 10; opts.rho_K = 10;  
fprintf(['ApriD with theta=',num2str(opts.theta),',alpha_K=',num2str(opts.alpha_K),',rho_K=',num2str(opts.rho_K),' is computing\n']); 
time = tic;
out_expect_ApriD_theta1_alpha10_rho100 = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD_theta1_alpha10_rho100 = toc(time);
% theta=10, calculated in Figure 4 
  
%%
savename = ['QCQP_exp_n_' num2str(n) '_p_' num2str(p) '_APriD_parameters'];
close all
clear time
save(savename) 

%% Plot Figures 4,5,6 in the paper with the above results.
plot_select_parameters_QCQP_expect
 