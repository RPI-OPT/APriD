function [out] = QCQP_expect_MSA(n,p,EHTH,EcTH,EcTc,EQ,Ea,Eb,X_min,X_max,opts)

if isfield(opts,'K')              K = opts.K;       else K = 1e5;      end
if isfield(opts,'alpha_K')  alpha_K = opts.alpha_K; else alpha_K = 1;  end
if isfield(opts,'rho_K')      rho_K = opts.rho_K;   else rho_K = 1;    end
 
if isfield(opts,'x')              x = opts.x;       else x =randn(1,n); end

if isfield(opts,'Jn')            Jn = opts.Jn;      else Jn = 10; end
if isfield(opts,'Jm')            Jm = opts.Jm;      else Jm = 10; end

if isfield(opts,'ks')         ks = opts.ks;      else    ks = (1:K+1);end

% def cal_grad_L_x;cal_grad_L_z ;

alpha = alpha_K/sqrt(K);
rho =  rho_K/sqrt(K);

z = 0;
%%% initial
xs = zeros(K+1,n);

num_ks = length(ks);
f0s = zeros(num_ks,1);
f1s = zeros(num_ks,1);
f0s_avgx = zeros(num_ks,1);
f1s_avgx = zeros(num_ks,1);
 
xs(1,:) = x; 
k_=0;
if ismember(1,ks)
    k_ = k_ + 1;
    f0s(k_) = cal_f0(x);
    f1s(k_) = cal_f1(x);
    f0s_avgx(k_) = f0s(k_);
    f1s_avgx(k_) = f1s(k_);
end

fprintf('QCQP_expect_MSA, Iteration:          ');
for k = 1+(1:K)
    
    fprintf('\b\b\b\b\b\b%6i',k);  
    
    u = cal_grad_L_x(x,z); %  use the new generated data
    w = cal_grad_L_z(x);    %  use the new generated data
    
    x = x - alpha*u.';
    x(x<X_min) = X_min;
    x(x>X_max) = X_max;
    z = max(z + rho*w,0);
    
    xs(k,:) = x;
    if ismember(k,ks)
        k_ = k_ + 1;
        f0s(k_) = cal_f0(x);
        f1s(k_) = cal_f1(x);
         
        x_average =  mean(xs(1:k,:),1);  
        f0s_avgx(k_) = cal_f0(x_average);
        f1s_avgx(k_) = cal_f1(x_average);
    end     
    
end
fprintf('\n');
%%
out.x = x;
out.f0s = f0s;
out.f1s = f1s;

out.x_average  =  x_average;
out.f0s_avgx = f0s_avgx;
out.f1s_avgx = f1s_avgx;
 
    function [Lx] = cal_grad_L_x(x,z)
        grad_f0 = zeros(n,1);
        for i = 1:Jn
            H_ = randn(p,n);  c_ = randn(p,1);  
            H_ = H_/norm(H_);c_ = c_/norm(c_);
            grad_f0 = grad_f0 + H_'*(H_*x'-c_);
        end
        grad_f0 = grad_f0/Jn;
        
        grad_f1 = 0;
        for i = 1:Jm
            Q_ = rand(n,n); Q_ = Q_*Q_'; Q_ = Q_/norm(Q_);
            a_ = randn(n,1);a_ = a_/norm(a_);
            
            grad_f1 = grad_f1 + (Q_*x'+a_);
        end
        grad_f1 = grad_f1/Jm;
        
        Lx =  grad_f0 + z*grad_f1;
    end

    function [Lz] =  cal_grad_L_z(x)
        f1 = 0;
        for i = 1:Jm
            Q_ = rand(n,n); Q_ = Q_*Q_'; Q_ = Q_/norm(Q_);
            a_ = randn(n,1);a_ = a_/norm(a_);
            b_ = 0.1 + rand();
            
            f1 = f1 + 1/2*x*Q_*x'+x*a_-b_;
        end
        Lz = f1/Jm;
    end
 
    function [f0] = cal_f0(x) 
        f0 = 0.5*(x*EHTH*x.'-2*EcTH*x.'+EcTc);
    end
    function [f1] = cal_f1(x)
        f1 = 0.5*x*EQ*x.'+x*Ea-Eb;
    end
end
