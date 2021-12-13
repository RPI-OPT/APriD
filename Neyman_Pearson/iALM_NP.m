function [w, out] = iALM_NP(X,y,alpha,opts)

% inexact Augmented Lagrangian method for Neyman-Pearson classification
% with logit loss function

[n,p] = size(X);

n_pos = sum(y==1);
n_neg = sum(y==-1);

X_pos = X(y==1,:);
X_neg = X(y==-1,:);
y_pos = y(y==1);
y_neg = y(y==-1);

u_pos = zeros(n_pos,1);
u_neg = zeros(n_neg,1);

if isfield(opts,'w0')       w0 = opts.w0;            else  w0 = randn(p,1); end
if isfield(opts,'maxit')    maxit = opts.maxit;      else  maxit = 100000;     end
if isfield(opts,'maxsubit') maxsubit = opts.maxsubit;else  maxsubit = 10000; end
if isfield(opts,'tol')      tol = opts.tol;          else  tol = 1e-10;      end
if isfield(opts,'beta')     beta = opts.beta;        else  beta = 1;        end


rho = beta;

inc = 1.5; dec = 2;

w = w0;
z = 0;
psi = 0;
eta0 = 1;
eta = eta0;

[grad_obj, gradPsi, f_obj, psi, f_constraint] = cal_grad_funVal(w,z);

obj0 = f_obj;
res0 = max(0, f_constraint);
hist_obj = f_obj;
hist_res = res0;
hist_obj_avg = f_obj;
hist_res_avg = res0;
hist_eta = eta0;
hist_z = z;

wavg = w;
sum_weight = 0;

k = 0;

num_grad = [];
num_f = [];

total_num_grad = 0;

%% main iterations
fprintf('iALM Iteration: (k, num_grad, total_num_grad)\n')
fprintf('                    ')
while total_num_grad < maxit
    k = k + 1;
    
    % update w
    
    %[w,subOut] = sub_solver(w0,z,beta,.001/k^2,maxsubit);
    [w,subOut] = sub_solver(w0,z,beta,max(.01/k^2,1e-8),maxsubit);
    
    num_grad = [num_grad; subOut.num_gradval];
    num_f = [num_f; subOut.num_fval];
    
    total_num_grad = total_num_grad + subOut.num_gradval;
    
    % update z
    
    z = z + rho*max(-z/beta, subOut.f_constraint);
    
    % save history
    obj = subOut.f_obj;
    res = max(0,subOut.f_constraint);
 
    
    hist_obj = [hist_obj; obj];
    hist_res = [hist_res; res];
    hist_eta = [hist_eta; eta];
    hist_z = [hist_z; z];    
    
    if norm(w-w0) < tol
        break;
    end
    
    w0 = w; 
    
    %fprintf('k = %i, num_grad = %i, total_num_grad = %i\n',k, num_grad(end),total_num_grad);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6i,%6i,%6i',k, num_grad(end),total_num_grad);
end
% end of main iteration

fprintf('\n\n');
out.hist_obj = hist_obj;
out.hist_res = hist_res;
out.hist_eta = hist_eta;
out.hist_z = hist_z;
out.num_grad = num_grad;
out.num_f = num_f;

%% sub_routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [w,out] = sub_solver(w0,z,beta,tol,submaxit)
        
        what = w0;
        
        w = w0;
        
        t0 = 1;
        num_fval = 0;
        num_gradval = 0;
        
        for subiter = 1:submaxit
            
            eta = eta0/dec;
            
            [grad_obj, gradPsi, f_obj, psi, ~] = cal_grad_funVal(what,z);
            
            gradL0 = gradPsi + grad_obj;
            Lval0 = psi + f_obj;
            
            num_fval = num_fval + 1;
            num_gradval = num_gradval + 1;
            
            Lval = inf;
            
            while Lval > Lval0 + gradL0'*(w-what) + .5*eta*norm(w-what)^2 + 100*eps
                
                eta = eta*inc;                
                w = what - gradL0/eta;
                [f_obj, psi] = cal_funVal(w,z);
                Lval = psi + f_obj;
                num_fval = num_fval + 1;                
            end
            
            
            [grad_obj, gradPsi, f_obj, psi, f_constraint] = cal_grad_funVal(w,z);
            
            if norm(grad_obj+gradPsi) < tol
                break;
            end
            
            t1 = (1+sqrt(1+4*t0^2))/2;
            what = w + (t0-1)/t1*(w-w0);
            t0 = t1;
            
            eta0 = eta; w0 = w;
            
        end
        
        out.num_fval = num_fval;
        out.num_gradval = num_gradval;
        out.f_obj = f_obj;
        out.f_constraint = f_constraint;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [f_obj, psi] = cal_funVal(w,z)
        
        expVal = -y_pos.*(X_pos*w);
        
        id1 = expVal>0;
        id2 = expVal<=0;
        
        expVal1 = expVal(id1);
        expVal2 = expVal(id2);
        
        u_pos(id1) = 1./(1+exp(-expVal1));
        u_pos(id2) = exp(expVal2)./(1+exp(expVal2));
        
        u_pos = -y_pos.*u_pos;
        
        f_obj = sum(expVal1)+sum(log(1+exp(-expVal1)))+sum(log(1+exp(expVal2)));
        f_obj = f_obj/n_pos;
        
        
        expVal = -y_neg.*(X_neg*w);
        
        id1 = expVal>0;
        id2 = expVal<=0;
        
        expVal1 = expVal(id1);
        expVal2 = expVal(id2);
        
        u_neg(id1) = 1./(1+exp(-expVal1));
        u_neg(id2) = exp(expVal2)./(1+exp(expVal2));
        
        u_neg = -y_neg.*u_neg;
        
        f_constraint = sum(expVal1)+sum(log(1+exp(-expVal1)))+sum(log(1+exp(expVal2)));
        f_constraint = f_constraint/n_neg - alpha;
        
        psi = -z^2/(2*beta);
        if beta*f_constraint + z > 0
            psi = z*f_constraint + .5*beta*f_constraint^2;
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [grad_obj, gradPsi, f_obj, psi, f_constraint] = cal_grad_funVal(w,z)
        
        expVal = -y_pos.*(X_pos*w);
        
        id1 = expVal>0;
        id2 = expVal<=0;
        
        expVal1 = expVal(id1);
        expVal2 = expVal(id2);
        
        u_pos(id1) = 1./(1+exp(-expVal1));
        u_pos(id2) = exp(expVal2)./(1+exp(expVal2));
        
        u_pos = -y_pos.*u_pos;
        grad_obj = X_pos'*u_pos/n_pos;
        
        f_obj = sum(expVal1)+sum(log(1+exp(-expVal1)))+sum(log(1+exp(expVal2)));
        f_obj = f_obj/n_pos;
        
        
        expVal = -y_neg.*(X_neg*w);
        
        id1 = expVal>0;
        id2 = expVal<=0;
        
        expVal1 = expVal(id1);
        expVal2 = expVal(id2);
        
        u_neg(id1) = 1./(1+exp(-expVal1));
        u_neg(id2) = exp(expVal2)./(1+exp(expVal2));
        
        u_neg = -y_neg.*u_neg;
        grad_constraint = X_neg'*u_neg/n_neg;
        
        f_constraint = sum(expVal1)+sum(log(1+exp(-expVal1)))+sum(log(1+exp(expVal2)));
        f_constraint = f_constraint/n_neg - alpha;
        
        psi = -z^2/(2*beta);
        gradPsi = zeros(p,1);
        if beta*f_constraint + z > 0
            psi = z*f_constraint + .5*beta*f_constraint^2;
            gradPsi = (beta*f_constraint + z)*grad_constraint;
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
