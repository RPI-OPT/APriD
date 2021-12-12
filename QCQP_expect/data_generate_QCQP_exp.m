function [EHTH,EcTH,EcTc,EQ,Ea,Eb] = data_generate_QCQP_exp(n,p,N,M)
Ea = zeros(n,1);
Eb = 0.6; 
EcTc = 1;

eHTH = zeros(n,n);  ecTH = zeros(1,n);  

fprintf('data generating, H:          ');
for i = 1:N
    if mod(i,100)==0
        fprintf('\b\b\b\b\b\b%6i',i);
    end
    H = randn(p,n);
    c = randn(p,1); 
    
    H = H/norm(H);
    c = c/norm(c);
    
    eHTH = eHTH + H.'*H;
    ecTH = ecTH + c.'*H;
end
EHTH =  eHTH/N;
EcTH =  ecTH/N; 
fprintf('\n');

fprintf('data generating, Q:          ');
EQ =  zeros(n,n);
for i = 1:M 
    if mod(i,100)==0
        fprintf('\b\b\b\b\b\b%6i',i);
    end
     Q = rand(n,n); Q = Q*Q';   Q = Q/norm(Q);  
     EQ = EQ+Q;
end
EQ = EQ/M;

fprintf('\n');
end