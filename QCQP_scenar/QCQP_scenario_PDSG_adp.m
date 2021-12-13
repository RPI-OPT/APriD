function [out] = QCQP_scenario_PDSG_adp(H,c,eHTH,ecTH,ecTc,Q,a,b,X_min,X_max,opts)

n = size(Q,1);

if isfield(opts,'K')              K = opts.K;       else K = 1e5;      end
if isfield(opts,'alpha_K')  alpha_K = opts.alpha_K; else alpha_K = 1;  end
if isfield(opts,'rho_K')      rho_K = opts.rho_K;   else rho_K = 1;    end 

if isfield(opts,'x')              x = opts.x;       else x =randn(1,n); end

if isfield(opts,'Jn')            Jn = opts.Jn;      else Jn = 10; end
if isfield(opts,'Jm')            Jm = opts.Jm;      else Jm = 10; end
 
if isfield(opts,'ks')         ks = opts.ks;         else ks = (1:K+1);end


if isfield(opts,'beta')        beta = opts.beta;    else beta = 1;   end 
if isfield(opts,'max_nrm')     max_nrm = opts.max_nrm;  else max_nrm = 1;   end
if isfield(opts,'one_eta2')    one_eta2 = opts.one_eta2; else one_eta2=10;    end
%% def f0; f1; 

%% step size
alpha = alpha_K/sqrt(K);
rho =  rho_K/sqrt(K);

z =zeros(1,length(b));
%% initial 

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

sum_grad = 0;

N_perm = randperm(size(H,3)); n0 = 1;
M_perm = randperm(size(Q,3)); m0 = 1;
fprintf('ada pdsg Iteration:          ');
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
  
%     u = cal_grad_L_x(x,z(m_selec),H_batch,c_batch,Q_batch,a_batch);
%     w = cal_grad_L_z(x,Q_batch,a_batch,b_batch);
    
    [u,w] = cal_grad_L_beta_xz(x,z(m_selec),beta,H_batch,c_batch,Q_batch,a_batch,b_batch);
    
    
%     x = x - alpha*u.';
    nrm_grad = norm(u);
    nrmlize_grad = u/max(max_nrm,nrm_grad);
    sum_grad = sum_grad + nrmlize_grad.^2; 
    x = x - (u./(sqrt(sum_grad/one_eta2)+1/alpha)).';
    
    x(x<X_min) = X_min;
    x(x>X_max) = X_max;
%     z(m_selec) = max(z(m_selec) + rho*w,0); 
    z(m_selec) = z(m_selec) + rho*max(-z(m_selec)/beta,w); 
    
    %%% keep record
    xs(k,:) = x; 
    if ismember(k,ks)
        k_ = k_ + 1;
        f0s(k_) = cal_f0(x);  
        f1s = cal_f1(x,Q,a,b);  
        f1s = max(f1s,0);
        f1s_avg(k_) = mean(f1s);
        f1s_max(k_) = max(f1s);
  
        x_average = mean(xs(1:k,:),1); 
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
 
% fprintf('ada pdsg, last x: f0, f1_avg, f1_max are')
% disp([f0s(end) f1s_avg(end) f1s_max(end)] )
% fprintf('ada pdsg, average x: f0, f1_avg, f1_max are')
% disp([f0s_avgx(end) f1s_avgx_avg(end) f1s_avgx_max(end)])
 
  
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

%     function [grad_f1] = cal_grad_f1(x,Q,a)  
%         grad_f1 = permute(sum(bsxfun(@times,Q,x),2),[1,3,2])+a;
%     end

%     function [grad_L_x] = cal_grad_L_x(x,z,H,c,Q,a) 
%         grad_L_x = cal_grad_f0(x,H,c) + mean(cal_grad_f1(x,Q,a).*z,2);
%     end

    function [u,w] = cal_grad_L_beta_xz(x,z,beta,H,c,Q,a,b)       
%         f1 =  cal_f1(x,Q,a,b) 
%         grad_f1 = cal_grad_f1(x,Q,a) 
        Qx = permute(sum(bsxfun(@times,Q,x),2),[1,3,2]);       
        f1 = 1/2*sum(bsxfun(@times,Qx,x'),1)+x*a-b';
        %grad_f1 = Qx+a; 
        
        u = cal_grad_f0(x,H,c) + mean( (Qx+a).*max(beta*f1+z,0) ,2);
        w = f1;
    end
   
end
