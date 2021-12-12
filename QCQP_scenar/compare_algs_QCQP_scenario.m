%% code for finite-sum structured QCQP with many constraints
%  need a given n before running

%%% min_{x in [X_min,X_max]^n} 1/(2N) sum_{i=1...N} \|H_i*x-c_i\|^2
%%% s.t.  1/2 x^TQ_jx+a_j^Tx-b_j<=0,j=1,...,M
%%
switch n
    case 10 
        p = 5;  
        opts.K = 50000; ks_max_intervial = 1000;
        fprintf(['finite-sum structured QCQP with many constraints with n=' num2str(n) ', p=' num2str(p) '\n']);
    case 200
        p = 150; 
        opts.K = 50000; ks_max_intervial = 100;
        fprintf(['finite-sum structured QCQP with many constraints with n=' num2str(n) ', p=' num2str(p) '\n']);
    otherwise
        error('new n, need give a new p');
end


N = 10000; M = 10000;
X_min = -10; X_max = 10;
% call for data generate
rng default
[H,c,Q,a,b,eHTH,ecTH,ecTc] = data_generate_QCQP_scenario(n,p,N,M); 

%%
opts.alpha_K = 10; opts.rho_K = sqrt(10); opts.eta = 0.04; opts.gam_K = 10;
opts.x = max(min(randn(1,n),X_max),X_min);

% record the results at the iteration  k  in  ks
ks= 1;
for k=1+(1:opts.K)
    if (mod(k,10^(floor(log10(k))))==0) ||  (mod(k,ks_max_intervial)==0) || k==opts.K+1
        ks=[ks,k];
    end
end
opts.ks=ks;

%%
fprintf('data generated; begin computing. \n')
%%% only work for n = 10; p = 5; for n = 200, p = 150; no results after 3hours
if n==10
    fprintf('begin cvx \n')
    time = tic;
    cvx_solver Mosek
    cvx_begin
    variable x(n)
    object =  0;
    
    minimize(1/2*(quad_form(x,eHTH)-2*ecTH*x+ecTc));
    for j = 1:M
        0.5*quad_form(x,Q(:,:,j))+a(:,j)'*x-b(j) <= 0;
    end
    X_min <= x <= X_max
    cvx_end
    
    time_scenario_cvx = toc(time);
    
    cvx_optval
    
    fprintf('end cvx. \n')
end

%%
fprintf('ada pdsg is computing\n'); 
time = tic;
opts.one_eta2=100; opts.alpha_K = 20; opts.rho_K =  sqrt(10); 
out_scenario_PDSG_adp = QCQP_scenario_PDSG_adp(H,c,eHTH,ecTH,ecTc,Q,a,b,X_min,X_max,opts);
time_scenario_PDSG_adp = toc(time);

%%
opts.alpha_K = 10; opts.rho_K = sqrt(10); opts.eta = 0.04; opts.gam_K = 10;
%%
fprintf('ApriD is computing\n'); 
time = tic;
out_scenario_ApriD = QCQP_scenario_ApriD(H,c,eHTH,ecTH,ecTc,Q,a,b,X_min,X_max,opts);
time_scenario_ApriD = toc(time);

%%
fprintf('MSA is computing\n'); 
time = tic;
out_scenario_MSA = QCQP_scenario_MSA(H,c,eHTH,ecTH,ecTc,Q,a,b,X_min,X_max,opts);
time_scenario_MSA = toc(time);
 
%%
fprintf('CSA is computing\n'); 
time = tic;
out_scenario_CSA = QCQP_scenario_CSA(H,c,eHTH,ecTH,ecTc,Q,a,b,X_min,X_max,opts);
time_scenario_CSA = toc(time);
 
%%
savename = ['QCQP_scenario_n_' num2str(n) '_p_' num2str(p) ...
    '_M_' num2str(M) '_N_' num2str(N) '_K_' num2str(opts.K) ...
    '_alpha_' num2str(opts.alpha_K) ...
    '_rho_' num2str(floor(opts.rho_K)) '_gam_' num2str(opts.gam_K)];

close all
clear H c Q a b time
save(savename)
%% fprintf(['QCQP scenario cvx takes time ' num2str(time_scenario_cvx) '\n']) 

fprintf(['computing time with n=' num2str(n) ', p=' num2str(p)  '\n'])
fprintf(['Example' ' & ' 'data or size'  ' & ' 'APriD' ' & ' 'MSA' ' & ' 'CSA' ' & '  'PDSG_adp'  '  \\\\\n'])
fprintf(['QCQP(5.4)' ' & ' '(' num2str(n) ',' num2str(p) ')' ' & ' num2str(time_scenario_ApriD,'%.1f') ' & ' num2str(time_scenario_MSA,'%.1f') ' & ' num2str(time_scenario_CSA,'%.1f') ' & ' num2str(time_scenario_PDSG_adp,'%.1f') '  \\\\\n'])

file_time = fopen('running_time.txt','a');
fprintf(file_time, ['QCQP(5.4)' ' & ' '(' num2str(n) ',' num2str(p) ')' ' & ' num2str(time_scenario_ApriD,'%.1f') ' & ' num2str(time_scenario_MSA,'%.1f') ' & ' num2str(time_scenario_CSA,'%.1f') ' & ' num2str(time_scenario_PDSG_adp,'%.1f') '  \\\\\n']);
fclose(file_time);

%% Plot Figure 3 in the paper with the above results.
plot_compare_algs_QCQP_scenario