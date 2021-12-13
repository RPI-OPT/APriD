close all
fprintf(['computing time on ' filename '\n'])
fprintf(['data set' ' & ' 'APriD' ' & ' 'MSA' ' & ' 'CSA' '  \\\\\n'])
fprintf([filename ' & ' num2str(time_ApriD,'%.1f') ' & ' num2str(time_MSA,'%.1f') ' & ' num2str(time_CSA,'%.1f') '  \\\\\n'])

time_max = max([time_ApriD,time_MSA,time_CSA ]);

%% 
figure('papersize',[20,4],'paperposition',[-1.7,0,23.2,4]);
set(gcf,'Position',[100 200 1600 300])
set(gca,'FontSize',20,'FontWeight','normal')
subplot(1,4,1)
semilogx(ks,0*ones(size(out_ApriD.f0s)),'k-', 'linewidth', 4);hold on 
semilogx(ks,abs(out_ApriD.f0s_avgx-f0_opt), 'b-', 'linewidth', 2);hold on
semilogx(ks,abs(out_MSA.f0s_avgx-f0_opt), 'g--', 'linewidth', 2);hold on
semilogx(ks,abs(out_CSA.f0s_avgx-f0_opt), 'r-.', 'linewidth', 2);hold on
semilogx(ks,abs(out_CSA.f0s_meanx-f0_opt), 'r-', 'linewidth', 2);hold off 
legend({'0','APriD','MSA','CSA1','CSA2'},'Location','southwest')
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20) 

subplot(1,4,3) 
semilogx(ks,f1_opt*ones(size(out_ApriD.f1s)),'k-', 'linewidth', 4);hold on
semilogx(ks,out_ApriD.f1s_avgx, 'b-', 'linewidth', 2);hold on
semilogx(ks,out_MSA.f1s_avgx, 'g--', 'linewidth', 2);hold on
semilogx(ks,out_CSA.f1s_avgx, 'r-.', 'linewidth', 2);hold on
semilogx(ks,out_CSA.f1s_meanx, 'r-', 'linewidth', 2);hold off
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('$f_1(\mathbf{x})$','interpreter','latex','FontSize',20)
 
subplot(1,4,2)
semilogx(time_max/ks(end)*ks,0*ones(size(out_ApriD.f0s)),'k-', 'linewidth', 4);hold on 
semilogx(time_CSA/ks(end)*ks,abs(out_ApriD.f0s_avgx-f0_opt), 'b-', 'linewidth', 2);hold on
semilogx(time_MSA/ks(end)*ks,abs(out_MSA.f0s_avgx-f0_opt), 'g--', 'linewidth', 2);hold on
semilogx(time_CSA/ks(end)*ks,abs(out_CSA.f0s_avgx-f0_opt), 'r-.', 'linewidth', 2);hold on
semilogx(time_CSA/ks(end)*ks,abs(out_CSA.f0s_meanx-f0_opt), 'r-', 'linewidth', 2);hold off  
xlabel('time','FontSize',20,'interpreter','latex') 
ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20) 
xlim([time_max/ks(end)*ks(1), time_max])

subplot(1,4,4) 
semilogx(time_max/ks(end)*ks,f1_opt*ones(size(out_ApriD.f1s)),'k-', 'linewidth', 4);hold on
semilogx(time_ApriD/ks(end)*ks,out_ApriD.f1s_avgx, 'b-', 'linewidth', 2);hold on
semilogx(time_MSA/ks(end)*ks,out_MSA.f1s_avgx, 'g--', 'linewidth', 2);hold on
semilogx(time_CSA/ks(end)*ks,out_CSA.f1s_avgx, 'r-.', 'linewidth', 2);hold on
semilogx(time_CSA/ks(end)*ks,out_CSA.f1s_meanx, 'r-', 'linewidth', 2);hold off
xlabel('time','FontSize',20,'interpreter','latex') 
ylabel('$f_1(\mathbf{x})$','interpreter','latex','FontSize',20)
xlim([time_max/ks(end)*ks(1), time_max])

sgtitle(filename,'FontSize',20)  
print(gcf,'-dpdf',savename);  
print(gcf,'-depsc',savename);  