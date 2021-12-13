%% code for QCQP in Expectation Form
%  need a given n before  running

%%% min_{x in [X_min,X_max]^n} E[1/2\|H_i*x-c_i\|^2]
%%% s.t.  E[1/2 x^TQ_jx+a_j^Tx-b_j]<=0, 

%%
switch n
    case 10 
        p=5;% n = 10; p = 5;
        opts.K = 50000;  ks_max_intervial = 10; 
        fprintf(['QCQP in Expectation Form with n=' num2str(n) ', p=' num2str(p) '\n']);
    case 200
        p=150;%n = 200; p = 150;
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
opts.alpha_K = 10; opts.rho_K = sqrt(10); opts.eta = 0.04; opts.gam_K = 10;


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
 
%%
fprintf('ApriD is computing\n'); 
time = tic;
out_expect_ApriD = QCQP_expect_ApriD(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_ApriD = toc(time);

%% 
fprintf('MSA is computing\n'); 
time = tic;
out_expect_MSA = QCQP_expect_MSA(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_MSA = toc(time); 
%%
fprintf('CSA is computing\n'); 
time = tic;
out_expect_CSA = QCQP_expect_CSA(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts);
time_expect_CSA = toc(time);

%%
savename = ['QCQP_exp_n_' num2str(n) '_p_' num2str(p) ...
    '_M_' num2str(M) '_N_' num2str(N) '_K_' num2str(opts.K) ...
    '_alpha_' num2str(opts.alpha_K) ...
    '_rho_' num2str(floor(opts.rho_K)) '_gam_' num2str(opts.gam_K)];
%
close all
clear time
save(savename)
%%
fprintf(['Example' ' & ' 'data or size' ' & ' 'APriD' ' & ' 'MSA' ' & ' 'CSA' ' & '  'PDSG_adp'  '  \\\\\n'])
fprintf(['QCQP(5.3)' ' & ' '(' num2str(n) ',' num2str(p) ')' ' & ' num2str(time_expect_ApriD,'%.1f') ' & ' num2str(time_expect_MSA,'%.1f') ' & ' num2str(time_expect_CSA,'%.1f') ' & ' ' - ' '  \\\\\n'])

file_time = fopen('running_time.txt','a');
fprintf(file_time, ['QCQP(5.3)' ' & ' '(' num2str(n) ',' num2str(p) ')' ' & ' num2str(time_expect_ApriD,'%.1f') ' & ' num2str(time_expect_MSA,'%.1f') ' & ' num2str(time_expect_CSA,'%.1f') ' & ' ' -'  '  \\\\\n']);
fclose(file_time);

%% Plot Figure 2 in the paper with the above results.
plot_compare_algs_QCQP_expect