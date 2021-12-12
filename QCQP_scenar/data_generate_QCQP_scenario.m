function [H,c,Q,a,b,eHTH,ecTH,ecTc] = data_generate_QCQP_scenario(n,p,N,M)

H = zeros(p,n,N);  eHTH = zeros(n,n);
c = zeros(p,N);    ecTH = zeros(1,n);  ecTc = 1;
fprintf('data generating, H:          ');
for i = 1:N
    if mod(i,100)==0
        fprintf('\b\b\b\b\b\b%6i',i);
    end
    H(:,:,i) = randn(p,n);
    c(:,i)   = randn(p,1);
     
    H(:,:,i) = H(:,:,i)/norm(H(:,:,i));
    c(:,i) = c(:,i)/norm(c(:,i));
    
    eHTH = eHTH + H(:,:,i).'*H(:,:,i);
    ecTH = ecTH + c(:,i).'*H(:,:,i); 
end
eHTH =  eHTH/N;
ecTH =  ecTH/N; 
fprintf('\n');

b = 0.1 + rand(M,1);
a = zeros(n,M);
Q = zeros(n,n,M);
fprintf('data generating, Q:          ');
for i = 1:M
    if mod(i,100)==0
        fprintf('\b\b\b\b\b\b%6i',i);
    end
    a(:,i) = randn(n,1);  a(:,i) = a(:,i)/norm(a(:,i));
    Q(:,:,i) = rand(n,n); Q(:,:,i) = Q(:,:,i)*Q(:,:,i)';Q(:,:,i) = Q(:,:,i)/norm(Q(:,:,i));
%     A = randn(n, n-5); Q(:,:,i) = A*A';Q(:,:,i) = Q(:,:,i)/norm(Q(:,:,i));
end
fprintf('\n');

end