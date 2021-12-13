function [out] =  QCQP_expect_CSA(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts)

if isfield(opts,'K')           K = opts.K;       else K = 1e5;      end

if isfield(opts,'eta')       eta = opts.eta;     else   eta = 1;    end
if isfield(opts,'gam_K')   gam_K = opts.gam_K;   else gam_K = 1;end

if isfield(opts,'x')           x = opts.x;       else x =randn(1,n); end

if isfield(opts,'Jg')         Jg = opts.Jg;      else Jg = 100;end
if isfield(opts,'Jn')         Jn = opts.Jn;      else Jn = 10; end
if isfield(opts,'Jm')         Jm = opts.Jm;      else Jm = 10; end

if isfield(opts,'ks')         ks = opts.ks;      else    ks = (1:K+1);end

% step size
gam =  gam_K/sqrt(K);
% eta

%%
% def f0; G; cal_grad_f0_x ;cal_grad_G_x ;
% and f1s which is used to record history information

%% initial
xs = zeros(K+1,n);

num_ks = length(ks);
f0s = zeros(num_ks,1);
f1s = zeros(num_ks,1);
f0s_avgx = zeros(num_ks,1);
f1s_avgx = zeros(num_ks,1);
f0s_meanx = zeros(num_ks,1);
f1s_meanx = zeros(num_ks,1);
 
xs(1,:) = x; 
k_=0;
if ismember(1,ks)
    k_ = k_ + 1;
    f0s(k_) = cal_f0(x);
    f1s(k_) = cal_f1(x);
    f0s_avgx(k_) = f0s(k_);
    f1s_avgx(k_) = f1s(k_); 
    f0s_meanx(1) = f0s(1);
    f1s_meanx(1) = f1s(1);
end   

tags = zeros(K+1,1); 
fprintf('QCQP_expect_CSA, Iteration:          ');
for k =  1+(1:K)
    fprintf('\b\b\b\b\b\b%6i',k);  
    
    g_estimated = cal_G(x);
     
    if g_estimated < eta
        h = cal_grad_f0(x);
        tags(k) = 1;
    else
        h = cal_grad_G(x);
    end
    
    x = x-gam*h.'; 
    x(x<X_min) = X_min;
    x(x>X_max) = X_max;  
    
    xs(k,:) = x;
    
    if ismember(k,ks)
        k_ = k_ + 1;
        f0s(k_) = cal_f0(x);
        f1s(k_) = cal_f1(x);
        
        x_average = mean(xs(logical(tags(1:k)),:),1);  
        f0s_avgx(k_) = cal_f0(x_average);
        f1s_avgx(k_) = cal_f1(x_average);
        
        x_mean = mean(xs(1:k,:),1);
        f0s_meanx(k_) = cal_f0(x_mean);
        f1s_meanx(k_) = cal_f1(x_mean);
    end 
    
end
fprintf('\n');
%%
out.x = x;
out.f0s = f0s;
out.f1s = f1s;

out.x_average = x_average;
out.f0s_avgx = f0s_avgx;
out.f1s_avgx = f1s_avgx; 

out.x_mean = x_mean;
out.f0s_meanx = f0s_meanx;
out.f1s_meanx = f1s_meanx; 

if sum(tags)==0
    out.average = ' "it alwasys update constraints!" ';
end

    function [grad_f0] = cal_grad_f0(x)
        grad_f0 = zeros(n,1);
        for i = 1:Jn
            H_ = randn(p,n);  c_ = randn(p,1);  
            H_ = H_/norm(H_);c_ = c_/norm(c_);
            grad_f0 = grad_f0 + H_'*(H_*x'-c_);
        end
        grad_f0 = grad_f0/Jn;
    end

    function [G] = cal_G(x)
        f1 = 0;
        for i = 1:Jg
            Q_ = rand(n,n); Q_ = Q_*Q_'; Q_ = Q_/norm(Q_);
            a_ = randn(n,1);a_ = a_/norm(a_);
            b_ = 0.1 + rand();
            
            f1 = f1 + 1/2*x*Q_*x'+x*a_-b_;
        end
        G = f1/Jg;
    end

    function [grad_G] = cal_grad_G(x)
        grad_f1 = 0;
        for i = 1:Jm
            Q_ = rand(n,n); Q_ = Q_*Q_'; Q_ = Q_/norm(Q_);
            a_ = randn(n,1);a_ = a_/norm(a_);
            
            grad_f1 = grad_f1 + (Q_*x'+a_);
        end
        grad_f1 = grad_f1/Jm;
        
        grad_G = grad_f1;
    end

    function [f0] = cal_f0(x) 
        f0 = 0.5*(x*EHTH*x.'-2*EcTH*x.'+EcTc);
    end
    function [f1] = cal_f1(x)
        f1 = 0.5*x*EQ*x.'+x*Ea-Eb;
    end

end
