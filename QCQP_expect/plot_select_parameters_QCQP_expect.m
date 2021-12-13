clc
close all
%% 
figure('papersize',[10,4],'paperposition',[-0.6,0,11.5,4]);
set(gcf,'Position',[100 200 800 300])
set(gca,'FontSize',20,'FontWeight','normal')
subplot(1,2,1) 

loglog(ks,abs(out_expect_ApriD10101000.f0s_avgx-cvx_optval),':','linewidth', 2,'DisplayName','$\rho$=100');hold on 
loglog(ks,abs(out_expect_ApriD1010100.f0s_avgx-cvx_optval),'-','linewidth', 2,'DisplayName','$\rho$=10');hold on 
loglog(ks,abs(out_expect_ApriD10101.f0s_avgx-cvx_optval),'-.','linewidth', 2,'DisplayName','$\rho$=1');hold on
loglog(ks,abs(out_expect_ApriD1010_1.f0s_avgx-cvx_optval),'--','linewidth', 2,'DisplayName','$\rho$=0.1');hold on
xlim([0,ks(end)])
legend('Location','southwest','interpreter','latex','FontSize',15)  
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20) 

subplot(1,2,2)
semilogx(ks,(out_expect_ApriD10101000.f1s_avgx), ':', 'linewidth', 2);hold on
semilogx(ks,(out_expect_ApriD1010100.f1s_avgx), '-', 'linewidth', 2);hold on 
semilogx(ks,(out_expect_ApriD10101.f1s_avgx), '-.', 'linewidth', 2);hold on
semilogx(ks,(out_expect_ApriD1010_1.f1s_avgx), '--', 'linewidth', 2);hold on
semilogx(ks,zeros(size(ks)), 'k-', 'linewidth', 1);hold off
xlim([0,ks(end)])
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$f_1(\mathbf{x})$','interpreter','latex','FontSize',20) 
 
sgtitle(['n = ' num2str(n) ', p = ' num2str(p)])
print(gcf,'-dpdf', ['QCQP_mosek_exp_n'  num2str(n) '_p'  num2str(p) '_fixed_theta_alpha']); 
print(gcf,'-depsc', ['QCQP_mosek_exp_n'  num2str(n) '_p'  num2str(p) '_fixed_theta_alpha']); 
 
%%
%%
figure('papersize',[10,4],'paperposition',[-0.6,0,11.5,4]);
set(gcf,'Position',[100 200 800 300])
subplot(1,2,1)  
semilogx(ks, (out_expect_ApriD10101.f0s_avgx-cvx_optval),':','linewidth', 2,'DisplayName','$\alpha=10,\rho=1$');hold on 
semilogx(ks, (out_expect_ApriD1052.f0s_avgx-cvx_optval),'-','linewidth', 2,'DisplayName','$\alpha=5,\rho=2$');hold on  
semilogx(ks, (out_expect_ApriD1025.f0s_avgx-cvx_optval),'-.','linewidth', 2,'DisplayName','$\alpha=2,\rho=5$');hold on 
semilogx(ks, (out_expect_ApriD10110.f0s_avgx-cvx_optval),'--','linewidth', 2,'DisplayName','$\alpha=1,\rho=10$'); hold off 
xlim([0,ks(end)])
legend('Location','northeast','interpreter','latex','FontSize',15)

xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20) 

subplot(1,2,2)
semilogx(ks,(out_expect_ApriD10101.f1s_avgx), ':', 'linewidth', 2);hold on
semilogx(ks,(out_expect_ApriD1052.f1s_avgx), '-', 'linewidth', 2);hold on
semilogx(ks,(out_expect_ApriD1025.f1s_avgx), '-.', 'linewidth', 2);hold on 
semilogx(ks,(out_expect_ApriD10110.f1s_avgx), '--', 'linewidth', 2);hold on
semilogx(ks,zeros(size(ks)), 'k-', 'linewidth', 1);hold off
xlim([0,ks(end)])
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$f_1(\mathbf{x})$','interpreter','latex','FontSize',20) 

sgtitle(['n = ' num2str(n) ', p = ' num2str(p)])
print(gcf,'-dpdf', ['QCQP_mosek_exp_n'  num2str(n) '_p'  num2str(p) '_fixed_theta']); 
print(gcf,'-depsc', ['QCQP_mosek_exp_n'  num2str(n) '_p'  num2str(p) '_fixed_theta']); 
 
%%
figure('papersize',[10,4],'paperposition',[-0.6,0,11.5,4]);
set(gcf,'Position',[100 200 800 300])
subplot(1,2,1)  
 
loglog(ks,abs(out_expect_ApriD0000110100.f0s_avgx-cvx_optval),'o-','MarkerSize',10,'MarkerIndices',[2,10,20,100,500,2000],'linewidth', 2,'DisplayName','$\theta=0.0001$');hold on  
loglog(ks,abs(out_expect_ApriD000110100.f0s_avgx-cvx_optval),'+-','MarkerSize',10,'MarkerIndices',2+[2,10,20,100,500,2000],'linewidth', 2,'DisplayName','$\theta=0.001$');hold on  
loglog(ks,abs(out_expect_ApriD00110100.f0s_avgx-cvx_optval),'*-','MarkerSize',10,'MarkerIndices',4+[2,10,20,100,500,2000],'linewidth', 2,'DisplayName','$\theta=0.01$');hold on  
loglog(ks,abs(out_expect_ApriD0110100.f0s_avgx-cvx_optval),'d-','MarkerSize',10,'MarkerIndices',6+[2,10,20,100,500,2000],'linewidth', 2,'DisplayName','$\theta=0.1$');hold on   
loglog(ks,abs(out_expect_ApriD110100.f0s_avgx-cvx_optval),'x-','MarkerSize',10,'MarkerIndices',7+[2,10,20,100,500,2000],'linewidth', 2,'DisplayName','$\theta=1$');hold on   
loglog(ks,abs(out_expect_ApriD1010100.f0s_avgx-cvx_optval),'s-','MarkerSize',10,'MarkerIndices',5+[2,10,20,100,500,2000],'linewidth', 2,'DisplayName','$\theta=10$');hold on  
xlim([0,ks(end)])
legend('Location','best','interpreter','latex','FontSize',15)
% title('$\alpha$=10,$\rho=10$','FontSize',20,'interpreter','latex')    
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20) 

subplot(1,2,2)
semilogx(ks,(out_expect_ApriD0000110100.f1s_avgx),'o-','MarkerSize',10,'MarkerIndices',[2,10,20,100,500,2000],'linewidth', 2);hold on
semilogx(ks,(out_expect_ApriD000110100.f1s_avgx), '+-','MarkerSize',10,'MarkerIndices',2+[2,10,20,100,500,2000], 'linewidth', 2);hold on
semilogx(ks,(out_expect_ApriD00110100.f1s_avgx), '*-','MarkerSize',10,'MarkerIndices',4+[2,10,20,100,500,2000], 'linewidth', 2);hold on
semilogx(ks,(out_expect_ApriD0110100.f1s_avgx), 'd-','MarkerSize',10,'MarkerIndices',6+[2,10,20,100,500,2000], 'linewidth', 2);hold on
semilogx(ks,(out_expect_ApriD110100.f1s_avgx), 'x-','MarkerSize',10,'MarkerIndices',7+[2,10,20,100,500,2000], 'linewidth', 2);hold on
semilogx(ks,(out_expect_ApriD1010100.f1s_avgx), 's-','MarkerSize',10,'MarkerIndices',5+[2,10,20,100,500,2000], 'linewidth', 2);hold on
semilogx(ks,zeros(size(ks)), 'k-', 'linewidth', 1);hold on 
xlim([0,ks(end)]) 
ylim([-0.7,1.7]) 
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$f_1(\mathbf{x})$','interpreter','latex','FontSize',20) 
 
sgtitle(['n = ' num2str(n) ', p = ' num2str(p)])
print(gcf,'-dpdf', ['QCQP_mosek_exp_n'  num2str(n) '_p'  num2str(p) '_vary_theta']); 
print(gcf,'-depsc', ['QCQP_mosek_exp_n'  num2str(n) '_p'  num2str(p) '_vary_theta']);

 