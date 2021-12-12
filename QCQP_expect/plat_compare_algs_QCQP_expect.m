close all
fprintf(['computing time with n=' num2str(n) ', p=' num2str(p)  '\n'])
fprintf(['(n,p)' ' & ' 'APriD' ' & ' 'MSA' ' & ' 'CSA' '  \\\\\n'])
fprintf(['(' num2str(n) ',' num2str(p) ')' ' & ' num2str(time_expect_ApriD,'%.1f') ' & ' num2str(time_expect_MSA,'%.1f') ' & ' num2str(time_expect_CSA,'%.1f') '  \\\\\n'])

time_max = max([time_expect_ApriD,time_expect_MSA,time_expect_CSA ]);

%%
% figure('papersize',[10,3],'paperposition',[0,0,10,3]);
figure('papersize',[20,4],'paperposition',[-1.7,0,23.2,4]);
set(gcf,'Position',[100 200 1600 300])
set(gca,'FontSize',20,'FontWeight','normal')
subplot(1,4,1) 
loglog(ks,abs(out_expect_ApriD.f0s_avgx-cvx_optval),'b-','linewidth', 2);hold on 
loglog(ks,abs(out_expect_MSA.f0s_avgx-cvx_optval),'g--','linewidth', 2);hold on 
loglog(ks,abs(out_expect_CSA.f0s_avgx-cvx_optval),'r-.','linewidth', 2);hold on 
loglog(ks,abs(out_expect_CSA.f0s_meanx-cvx_optval),'r-','linewidth', 2);hold off 
legend({'APriD','MSA','CSA1','CSA2'},'Location','southwest')
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20) 
% xlim([ks(1), 1e5])

subplot(1,4,3)
semilogx(zeros(1,5e4), 'k-', 'linewidth', 1);hold on
semilogx(ks,(out_expect_ApriD.f1s_avgx), 'b-', 'linewidth', 2);hold on
semilogx(ks,(out_expect_MSA.f1s_avgx), 'g--', 'linewidth', 2);hold on
semilogx(ks,(out_expect_CSA.f1s_avgx), 'r-.', 'linewidth', 2);hold on
semilogx(ks,(out_expect_CSA.f1s_meanx), 'r-', 'linewidth', 2);hold off
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$f_1(\mathbf{x})$','interpreter','latex','FontSize',20) 
% xlim([ks(1), ks(end)])
% figure('papersize',[10,3],'paperposition',[0,0,10,3]);
% figure('papersize',[10,4],'paperposition',[-0.6,0,11.5,4]);
% set(gcf,'Position',[100 200 800 300])
% set(gca,'FontSize',20,'FontWeight','normal')
subplot(1,4,2)
loglog(time_expect_ApriD/ks(end)*ks,abs(out_expect_ApriD.f0s_avgx-cvx_optval),'b-','linewidth', 2);hold on 
loglog(time_expect_MSA/ks(end)*ks,abs(out_expect_MSA.f0s_avgx-cvx_optval),'g--','linewidth', 2);hold on 
loglog(time_expect_CSA/ks(end)*ks,abs(out_expect_CSA.f0s_avgx-cvx_optval),'r-.','linewidth', 2);hold on 
loglog(time_expect_CSA/ks(end)*ks,abs(out_expect_CSA.f0s_meanx-cvx_optval),'r-','linewidth', 2);hold off 
% legend({'APriD','MSA','CSA1','CSA2'},'Location','southwest') 
xlabel('time','FontSize',20,'interpreter','latex') 
ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20) 
xlim([time_max/ks(end)*ks(1), time_max])

subplot(1,4,4) 
semilogx(time_max/ks(end)*ks, zeros(1,size(ks,2)), 'k-', 'linewidth', 1);hold on
semilogx(time_expect_ApriD/ks(end)*ks,(out_expect_ApriD.f1s_avgx), 'b-', 'linewidth', 2);hold on
semilogx(time_expect_MSA/ks(end)*ks,(out_expect_MSA.f1s_avgx), 'g--', 'linewidth', 2);hold on
semilogx(time_expect_CSA/ks(end)*ks,(out_expect_CSA.f1s_avgx), 'r-.', 'linewidth', 2);hold on
semilogx(time_expect_CSA/ks(end)*ks,(out_expect_CSA.f1s_meanx), 'r-', 'linewidth', 2);hold off
xlabel('time','FontSize',20,'interpreter','latex') 
ylabel('$f_1(\mathbf{x})$','interpreter','latex','FontSize',20)
xlim([time_max/ks(end)*ks(1), time_max])

sgtitle(['n = ' num2str(n) ', p = ' num2str(p)])

print(gcf,'-dpdf',savename); 
print(gcf,'-depsc', savename);