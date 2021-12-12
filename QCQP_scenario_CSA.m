function [out] = QCQP_scenario_CSA(H,c,eHTH,ecTH,ecTc,Q,a,b,X_min,X_max,opts)

n = size(Q,1);

if isfield(opts,'K')           K = opts.K;       else K = 1e5;      end

if isfield(opts,'eta')       eta = opts.eta;     else   eta = 1;    end
if isfield(opts,'gam_K')   gam_K = opts.gam_K;   else gam_K = 1/200;end

if isfield(opts,'x')           x = opts.x;       else x =randn(1,n); end

if isfield(opts,'Jg')         Jg = opts.Jg;      else Jg = 100;end
if isfield(opts,'Jn')         Jn = opts.Jn;      else Jn = 10; end
if isfield(opts,'Jm')         Jm = opts.Jm;      else Jm = 10; end

if isfield(opts,'ks')         ks = opts.ks;      else    ks = (1:K+1);end

%% def f0; G; cal_grad_f0_x ;cal_grad_G_x ;
% and  f1s which is used to record history information

%% step size
gam =  gam_K/sqrt(K);
% eta

%% initial
xs = zeros(K+1,n);
num_ks = length(ks);  
f0s = zeros(num_ks,1);
f1s_avg = zeros(num_ks,1);
f1s_max = zeros(num_ks,1);
f0s_avgx = zeros(num_ks,1);
f1s_avgx_avg = zeros(num_ks,1);
f1s_avgx_max = zeros(num_ks,1);
f0s_meanx = zeros(num_ks,1);
f1s_meanx_avg = zeros(num_ks,1);
f1s_meanx_max = zeros(num_ks,1);

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
    f0s_meanx(k_) = f0s(1);
    f1s_meanx_avg(k_) = f1s_avg(1);
    f1s_meanx_max(k_) = f1s_max(1);
end 

tags = zeros(K+1,1); 
 
g_perm = randperm(size(Q,3)); g0 = 1;
N_perm = randperm(size(H,3)); n0 = 1;
M_perm = randperm(size(Q,3)); m0 = 1;
fprintf('CSA  Iteration:          ');
for k = 1+(1:K)
    fprintf('\b\b\b\b\b\b%6i',k);
    % randomly select Jg data from constrain
    g1 = g0 + Jg - 1;
    if g1 > size(Q,3)
        g_perm = randperm(size(Q,3)); g0 = 1;
        g1 = g0 + Jg - 1;
    end
    g_selec = g_perm(g0:g1); g0 = g1+1;
    
    g_estimated = cal_G(x,Q(:,:,g_selec),a(:,g_selec),b(g_selec));
    
    % cyclic select x
    if g_estimated < eta
        n1 = n0 + Jn - 1;
        if n1 > size(H,3)
            N_perm = randperm(size(H,3)); n0 = 1;
            n1 = n0 + Jn - 1;
        end
        n_selec = N_perm(n0:n1); n0 = n1+1;
        
        h = cal_grad_f0(x,H(:,:,n_selec),c(:,n_selec));
        
        tags(k) = 1;
    else
        m1 = m0 + Jm - 1;
        if m1 > size(H,3)
            M_perm = randperm(size(Q,3)); m0 = 1;
            m1 = m0 + Jm - 1;
        end
        m_selec = M_perm(m0:m1); m0 = m1+1;
        
        h = cal_grad_G(x,Q(:,:,m_selec),a(:,m_selec),b(m_selec));
    end
    
    x = x-gam*h.';
    
    x(x<X_min) = X_min;
    x(x>X_max) = X_max; 
    
    %%% keep record
    xs(k,:) = x;
    if ismember(k,ks)
         k_ = k_ + 1;
         f0s(k_) = cal_f0(x); 
         f1s = cal_f1(x,Q,a,b);
         f1s = max(f1s,0);
         f1s_avg(k_) = mean(f1s);
         f1s_max(k_) = max(f1s);
         
         x_average = mean(xs(logical(tags(1:k)),:),1);
         % at the begining, there would be NaN
         f0s_avgx(k_) = cal_f0(x_average);  
         f1s_avgx = max(cal_f1(x_average,Q,a,b),0);
         f1s_avgx_avg(k_) = mean(f1s_avgx);
         f1s_avgx_max(k_) = max(f1s_avgx);
         
         x_mean = mean(xs(1:k,:),1);
         f0s_meanx(k_) = cal_f0(x_mean);%cal_f0(x_average,H,c);
         f1s_meanx = max(cal_f1(x_mean,Q,a,b),0);
         f1s_meanx_avg(k_) = mean(f1s_meanx);
         f1s_meanx_max(k_) = max(f1s_meanx);
    end
end
fprintf('\n');
%%
out.x = x;
out.f0s = f0s;
out.f1s_avg = f1s_avg;
out.f1s_max = f1s_max;

out.f0s_avgx = f0s_avgx;
out.f1s_avgx_avg = f1s_avgx_avg;
out.f1s_avgx_max = f1s_avgx_max;

out.f0s_meanx = f0s_meanx;
out.f1s_meanx_avg = f1s_meanx_avg;
out.f1s_meanx_max = f1s_meanx_max;
 
if sum(tags)==0
    out.average = ' "it alwasys update constraints!" ';
end
  
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

    function [G] = cal_G(x,Q,a,b)
        f1 = cal_f1(x,Q,a,b); 
        G = mean(f1.*(f1>=0));
    end
   
    function [grad_f1] = cal_grad_f1(x,Q,a) 
        grad_f1 = permute(sum(bsxfun(@times,Q,x),2),[1,3,2])+a;
    end

    function [grad_G] = cal_grad_G(x,Q,a,b) 
        f1 = cal_f1(x,Q,a,b); 
        grad_f1 =  cal_grad_f1(x,Q,a);
        grad_G  = mean(grad_f1.*(f1>=0),2);  
    end

end
