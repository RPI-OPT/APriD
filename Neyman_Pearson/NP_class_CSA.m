function [out] = NP_class_CSA(M_p,M_n,opts)

[m_p, d]= size(M_p);
m_n  = size(M_n,1);

if isfield(opts,'K')           K = opts.K;       else     K = 1e5;  end
if isfield(opts,'eta')       eta = opts.eta;     else   eta = 1;    end
if isfield(opts,'gam_K')   gam_K = opts.gam_K;   else gam_K = 1/200;end
if isfield(opts,'f1_r')     f1_r = opts.f1_r;    else  f1_r = 0.4;  end
if isfield(opts,'x')           x = opts.x;       else     x =randn(1,d);end
if isfield(opts,'Jg')         Jg = opts.Jg;      else    Jg = min(floor(m_n/10),100);end
if isfield(opts,'Jp')         Jp = opts.Jp;      else    Jp = 10;end
if isfield(opts,'Jn')         Jn = opts.Jn;      else    Jn = 10;end

if isfield(opts,'ks')         ks = opts.ks;      else    ks = (1:K+1);end
num_ks = length(ks);

gam =  gam_K/sqrt(K);
%eta = 4*gam;

f0 = @(x)     mean(log(one_add_e_z(-M_p*x')));
f1 = @(M_n,x) mean(log(one_add_e_z( M_n*x'))) - f1_r;
 
grad_f0 = @(x_p,x) -mean(bsxfun(@rdivide,x_p,one_add_e_z(x_p*x'))); 
grad_f1 = @(x_n,x) mean(bsxfun(@rdivide, x_n, one_add_e_z(-x_n*x')));

%%% initial
xs = zeros(K+1,d);
tags = zeros(K+1,1);

f0s = zeros(num_ks,1);
f1s = zeros(num_ks,1);
f0s_avgx = zeros(num_ks,1);
f1s_avgx = zeros(num_ks,1);
f0s_meanx = zeros(num_ks,1);
f1s_meanx = zeros(num_ks,1);

xs(1,:)=x;
k_=0;
if ismember(1,ks)
     k_ = k_ + 1;
    f0s(k_) = f0(x);
    f1s(k_) = f1(M_n,x);
    f0s_avgx(k_) = f0s(k_);
    f1s_avgx(k_) = f1s(k_);
    f0s_meanx(k_) = f0s(k_);
    f1s_meanx(k_) = f1s(k_);
end

g_perm = randperm(m_n); g0 = 1;
p_perm = randperm(m_p); p0 = 1;
n_perm = randperm(m_n); n0 = 1;
fprintf('CSA Iteration:          ');
for k = 1+(1:K)
      
    fprintf('\b\b\b\b\b\b%6i',k);
    g1 = g0 + Jg - 1;
    if g1 > m_n
        g_perm = randperm(m_n); g0 = 1;
        g1 = g0 + Jg - 1;
    end
    g_selec = g_perm(g0:g1); g0 = g1+1;
    
    g_estimated = f1(M_n(g_selec,:),x); 
    
    if g_estimated < eta
        p1 = p0 + Jp - 1;
        if p1 > m_p
            p_perm = randperm(m_p); p0 = 1;
            p1 = p0 + Jp - 1;
        end
        p_selec = p_perm(p0:p1); p0 = p1+1;
        
        x_p = M_p(p_selec,:);% positive sample
        
        h = grad_f0(x_p,x);
        tags(k) = 1;
    else
        n1 = n0 + Jn - 1;
        if n1 > m_n
            n_perm = randperm(m_n); n0 = 1;
            n1 = n0 + Jn - 1;
        end
        n_selec = n_perm(n0:n1); n0 = n1+1;
        
        x_n = M_n(n_selec,:);% negative sample
        
        h = grad_f1(x_n,x);
    end
    
    x = x-gam*h;
    
    xs(k,:) = x;
    
    if ismember(k,ks)
        k_ = k_ + 1;
        f0s(k_) = f0(x);
        f1s(k_) = f1(M_n,x);
        
        x_average = mean(xs(logical(tags(1:k)),:),1);
        f0s_avgx(k_) = f0(x_average);
        f1s_avgx(k_) = f1(M_n,x_average);
        
        x_mean = mean(xs(1:k,:),1);
        f0s_meanx(k_) = f0(x_mean);
        f1s_meanx(k_) = f1(M_n,x_mean);
    end
end
fprintf('\n');

out.x = x;
out.f0s = f0s;
out.f1s = f1s;

acc_p = mean(M_p*x'>0);
acc_n = mean(M_n*x'<0);
acc_all = mean([M_p*x'>0;M_n*x'<0]);
out.last = [acc_all acc_p acc_n];

out.x_average = x_average;
out.f0s_avgx = f0s_avgx;
out.f1s_avgx = f1s_avgx; 

if sum(logical(tags))>=1   
    acc_p = mean(M_p*x_average'>0);
    acc_n = mean(M_n*x_average'<0);
    acc_all = mean([M_p*x_average'>0;M_n*x_average'<0]);
    out.average = [acc_all acc_p acc_n];
else
    out.average = ' "it alwasys update constraints!" ';
end

out.x_mean = x_mean;
out.f0s_meanx = f0s_meanx;
out.f1s_meanx = f1s_meanx; 
 
acc_p = mean(M_p*x_mean'>0);
acc_n = mean(M_n*x_mean'<0);
acc_all = mean([M_p*x_mean'>0;M_n*x_mean'<0]);
out.mean = [acc_all acc_p acc_n];
 
    function ez_1 = one_add_e_z(z)
        ez_1 = zeros(size(z));
        ez_1(z<=0) = 1+exp(z(z<=0));
        ez_1(z>0) = (1+exp(-z(z>0)))./exp(-z(z>0)); 
    end
end
