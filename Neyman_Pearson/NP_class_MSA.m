function [out] = NP_class_MSA(M_p,M_n,opts)

d = size(M_p,2);

if isfield(opts,'K')              K = opts.K;       else K = 1e5;      end
if isfield(opts,'alpha_K')  alpha_K = opts.alpha_K; else alpha_K = 1;  end
if isfield(opts,'rho_K')      rho_K = opts.rho_K;   else rho_K = 1;    end 

if isfield(opts,'f1_r')        f1_r = opts.f1_r;    else f1_r = 0.4;   end
if isfield(opts,'x')              x = opts.x;       else x =randn(1,d);end
if isfield(opts,'z')              z = opts.z;       else z = 0;        end

if isfield(opts,'Jp')         Jp = opts.Jp;      else    Jp = 10;end
if isfield(opts,'Jn')         Jn = opts.Jn;      else    Jn = 10;end

if isfield(opts,'ks')         ks = opts.ks;      else    ks = (1:K+1);end
num_ks = length(ks);


alpha = alpha_K/sqrt(K); 
rho =  rho_K/sqrt(K);

m_p = size(M_p,1);
m_n = size(M_n,1);

f0 =@(x) mean(log(one_add_e_z(-x*M_p')));
f1 =@(x) mean(log(one_add_e_z( x*M_n'))) - f1_r; 
 
grad_f0 = @(x_p,x) -mean(bsxfun(@rdivide,x_p,one_add_e_z(x_p*x')));  
grad_f1 = @(x_n,x) mean(bsxfun(@rdivide, x_n, one_add_e_z(-x_n*x')));
grad_L_x = @(x_p,x_n,x,z) grad_f0(x_p,x) + z*grad_f1(x_n,x);

%grad_L_z = @(x_n,x) mean(log(one_add_e_z(x*x_n')))-f1_r;
grad_L_z = @(x_n,x) mean(log(one_add_e_z(x_n*x'))) - f1_r;

%%% initial
xs = zeros(K+1,d);
f0s = zeros(num_ks,1);
f1s = zeros(num_ks,1);
f0s_avgx = zeros(num_ks,1);
f1s_avgx = zeros(num_ks,1);

xs(1,:)=x;
k_=0;
if ismember(1,ks)
     k_ = k_ + 1;
    f0s(k_) = f0(x);
    f1s(k_) = f1(x);
    f0s_avgx(k_) = f0s(k_);
    f1s_avgx(k_) = f1s(k_);
end

p_perm = randperm(m_p); p0 = 1;
n_perm = randperm(m_n); n0 = 1;
fprintf('MSA Iteration:          ');
for k = 1+(1:K)
    
    fprintf('\b\b\b\b\b\b%6i',k);
    
    p1 = p0 + Jp - 1;
    if p1 > m_p
        p_perm = randperm(m_p); p0 = 1;
        p1 = p0 + Jp - 1;
    end
    p_selec = p_perm(p0:p1); p0 = p1+1;
    
    n1 = n0 + Jn - 1;
    if n1 > m_n
        n_perm = randperm(m_n); n0 = 1;
        n1 = n0 + Jn - 1;
    end
    n_selec = n_perm(n0:n1); n0 = n1+1; 
    
    x_p = M_p(p_selec,:);% positive sample
    x_n = M_n(n_selec,:);% negative sample 
    
    u = grad_L_x(x_p,x_n,x,z);
    w = grad_L_z(x_n,x); 
    
    x = x - alpha*u;
    z = max(z + rho*w,0);
    
    xs(k,:) = x;  
    if ismember(k,ks)
        k_ = k_ + 1;
        f0s(k_) = f0(x);
        f1s(k_) = f1(x);
        
        x_average = mean(xs(1:k,:),1);
        f0s_avgx(k_) = f0(x_average);
        f1s_avgx(k_) = f1(x_average);
    end
 
end
fprintf('\n'); 
%%
out.x = x;
out.f0s = f0s;
out.f1s = f1s;

acc_p = mean(M_p*x'>0);
acc_n = mean(M_n*x'<0);
acc_all = mean([M_p*x'>0;M_n*x'<0]);
out.last = [acc_all acc_p acc_n];

out.x_average  =  x_average;
out.f0s_avgx = f0s_avgx;
out.f1s_avgx = f1s_avgx;
 
acc_p = mean(M_p*x_average'>0);
acc_n = mean(M_n*x_average'<0);
acc_all = mean([M_p*x_average'>0;M_n*x_average'<0]);
out.average = [acc_all acc_p acc_n];
 
    function ez_1 = one_add_e_z(z)
        ez_1 = zeros(size(z));
        ez_1(z<=0) = 1+exp(z(z<=0));
        ez_1(z>0) = (1+exp(-z(z>0)))./exp(-z(z>0)); 
    end
end
