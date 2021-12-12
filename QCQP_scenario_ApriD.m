function [out] = QCQP_scenario_ApriD(H,c,eHTH,ecTH,ecTc,Q,a,b,X_min,X_max,opts)

n = size(Q,1);

if isfield(opts,'K')              K = opts.K;       else K = 1e5;      end
if isfield(opts,'alpha_K')  alpha_K = opts.alpha_K; else alpha_K = 1;  end
if isfield(opts,'rho_K')      rho_K = opts.rho_K;   else rho_K = 1;    end

if isfield(opts,'beta1')      beta1 = opts.beta1;   else beta1 = 0.9;  end
if isfield(opts,'beta2')      beta2 = opts.beta2;   else beta2 = 0.99; end

if isfield(opts,'x')              x = opts.x;       else x =randn(1,n); end

if isfield(opts,'Jn')            Jn = opts.Jn;      else Jn = 10; end
if isfield(opts,'Jm')            Jm = opts.Jm;      else Jm = 10; end

if isfield(opts,'ks')         ks = opts.ks;      else    ks = (1:K+1);end
 
theta = 10;

%% def f0; f1; cal_grad_L_x ;
cal_grad_L_z = @(x,Q,a,b) cal_f1(x,Q,a,b);

%% step size
alpha = alpha_K/sqrt(K);
rho =  rho_K/sqrt(K);

z =zeros(1,length(b));
%% initial
m = zeros(size(x));
v = zeros(size(x));
v_hat = zeros(size(x));

xs = zeros(K+1,n);

num_ks = length(ks);  
f0s = zeros(num_ks,1);
f1s_avg = zeros(num_ks,1);
f1s_max = zeros(num_ks,1);
f0s_avgx = zeros(num_ks,1);
f1s_avgx_avg = zeros(num_ks,1);
f1s_avgx_max = zeros(num_ks,1);

xs(1,:) = x; 
k_=0;
if ismember(1,ks)
    k_ = k_ + 1;
    f0s(k_) = cal_f0(x);  
    f1s = max(cal_f1(x,Q,a,b),0); 
    
    f1s_avg(k_) = mean(f1s);
    f1s_max(k_) = max(f1s);
    f0s_avgx(k_) = f0s(1);
    f1s_avgx_avg(k_) = f1s_avg(1);
    f1s_avgx_max(k_) = f1s_max(1);
end
 
N_perm = randperm(size(H,3)); n0 = 1;
M_perm = randperm(size(Q,3)); m0 = 1;
fprintf('ApriD Iteration:          ');
for k = 1+(1:K)
    
    fprintf('\b\b\b\b\b\b%6i',k);
    
    n1 = n0 + Jn - 1;
    if n1 > size(H,3)
        N_perm = randperm(size(H,3)); n0 = 1;
        n1 = n0 + Jn - 1;
    end
    n_selec = N_perm(n0:n1); n0 = n1+1;
    
    m1 = m0 + Jm - 1;
    if m1 > size(H,3)
        M_perm = randperm(size(Q,3)); m0 = 1;
        m1 = m0 + Jm - 1;
    end
    m_selec = M_perm(m0:m1); m0 = m1+1;
    
    H_batch = H(:,:,n_selec);  c_batch = c(:,n_selec);  
    Q_batch = Q(:,:,m_selec);  a_batch = a(:,m_selec);  b_batch = b(m_selec);  
  
    u = cal_grad_L_x(x,z(m_selec),H_batch,c_batch,Q_batch,a_batch);
    w = cal_grad_L_z(x,Q_batch,a_batch,b_batch);
    
    m = beta1*m+(1-beta1)*u.';
    v = beta2*v+(1-beta2)*(u.').^2/max(1,norm(u)^2/theta);
    v_hat = max(v,v_hat);
    
    x = x - alpha*m./(sqrt(v_hat)+1e-15);
    x(x<X_min) = X_min;
    x(x>X_max) = X_max;
    z(m_selec) = max(z(m_selec) + rho*w,0);
    
    %%% keep record
    xs(k,:) = x;
    if ismember(k,ks)
        k_ = k_ + 1;
        f0s(k_) = cal_f0(x);  
        f1s = cal_f1(x,Q,a,b);  
        f1s = max(f1s,0);
        f1s_avg(k_) = mean(f1s);
        f1s_max(k_) = max(f1s);
  
        beta1_hat = 1 - beta1.^((k):-1:1);
        x_average = beta1_hat/sum(beta1_hat)*xs(1:k,:); 
        f0s_avgx(k_) = cal_f0(x_average); 
        f1s_avgx = cal_f1(x_average,Q,a,b); 
        f1s_avgx = max(f1s_avgx,0);
        f1s_avgx_avg(k_) = mean(f1s_avgx);
        f1s_avgx_max(k_) = max(f1s_avgx);
    end
end
fprintf('\n');
%%
out.x = x;
out.x_average =  x_average;
out.f0s = f0s;
out.f1s_avg = f1s_avg;
out.f1s_max = f1s_max;

out.f0s_avgx = f0s_avgx;
out.f1s_avgx_avg = f1s_avgx_avg;
out.f1s_avgx_max = f1s_avgx_max;
  
    function [f0] =cal_f0(x) 
        f0  = 1/2*(x*eHTH*x'-2*ecTH*x'+ecTc);
    end

    function [grad_f0] = cal_grad_f0(x,H,c)
        [p,N] = size(c);
        Hx_c = permute(sum(bsxfun(@times, H,x),2),[1,3,2])-c;
        HHx_c = sum(bsxfun(@times,H, reshape(Hx_c,[p,1,N])),1);
        grad_f0 = mean(HHx_c,3)';
    end

    function [f1] = cal_f1(x,Q,a,b)  
        f1 = 1/2*sum(bsxfun(@times,permute(sum(bsxfun(@times,Q,x),2),[1,3,2]),x'),1)+x*a-b';
    end

    function [grad_f1] = cal_grad_f1(x,Q,a) 
        grad_f1 = permute(sum(bsxfun(@times,Q,x),2),[1,3,2])+a;
    end

    function [grad_L_x] = cal_grad_L_x(x,z,H,c,Q,a) 
        grad_L_x = cal_grad_f0(x,H,c) + mean(cal_grad_f1(x,Q,a).*z,2);
    end 

end
