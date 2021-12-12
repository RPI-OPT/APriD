close all
fprintf(['computing time with n=' num2str(n) ', p=' num2str(p)  '\n'])
fprintf(['(n,p)' ' & ' 'APriD' ' & ' 'MSA' ' & ' 'CSA' ' & '  'PDSG_adp'  '  \\\\\n'])
fprintf(['(' num2str(n) ',' num2str(p) ')' ' & ' num2str(time_scenario_ApriD,'%.1f') ' & ' num2str(time_scenario_MSA,'%.1f') ' & ' num2str(time_scenario_CSA,'%.1f') ' & ' num2str(time_scenario_ada_pdsg,'%.1f') '  \\\\\n'])

time_max=max([time_scenario_ApriD,time_scenario_MSA,time_scenario_CSA,time_scenario_ada_pdsg ]);
%%
figure('papersize',[15,8],'paperposition',[0,0,15,8]);
set(gcf,'Position',[100 200 1200 600])
subplot(2,3,1)
if n==10
    loglog(ks,abs(out_scenario_ApriD.f0s_avgx-cvx_optval), 'b-', 'linewidth', 2);hold on
    loglog(ks,abs(out_scenario_MSA.f0s_avgx-cvx_optval), 'g--', 'linewidth', 2);hold on
    loglog(ks,abs(out_scenario_CSA.f0s_avgx-cvx_optval), 'r-.', 'linewidth', 2);hold on
    loglog(ks,abs(out_scenario_CSA.f0s_meanx-cvx_optval), 'r-', 'linewidth', 2);hold on 
    loglog(ks,abs(out_scenario_ada_pdsg.f0s_avgx-cvx_optval), 'k-', 'linewidth', 2);hold off 
    ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20)
end
if n==200
    optval_estimate = min([
        min(out_scenario_ApriD.f0s_avgx(out_scenario_ApriD.f1s_avgx_max<=0)),...
        min(out_scenario_MSA.f0s_avgx(out_scenario_MSA.f1s_avgx_max<=0)),...
        min(out_scenario_CSA.f0s_avgx(out_scenario_CSA.f1s_avgx_max<=0)),... 
        min(out_scenario_ada_pdsg.f0s_avgx(out_scenario_ada_pdsg.f1s_avgx_max<=0))...
        ]);
    loglog(ks,out_scenario_ApriD.f0s_avgx-optval_estimate, 'b-', 'linewidth', 2);hold on
    loglog(ks,out_scenario_MSA.f0s_avgx-optval_estimate, 'g--', 'linewidth', 2);hold on
    loglog(ks,out_scenario_CSA.f0s_avgx-optval_estimate, 'r-.', 'linewidth', 2);hold on
    loglog(ks,out_scenario_CSA.f0s_meanx-optval_estimate, 'r-', 'linewidth', 2);hold on
%     loglog(ks,out_scenario_pdsg.f0s_avgx-optval_estimate, 'k--', 'linewidth', 2);hold on
    loglog(ks,out_scenario_ada_pdsg.f0s_avgx-optval_estimate, 'k-', 'linewidth', 2);hold on
%     loglog(ks,out_scenario_ada_pdsg1.f0s_avgx-optval_estimate, 'k-.', 'linewidth', 2);hold on
%     loglog(ks,out_scenario_ada_pdsg2.f0s_avgx-optval_estimate, 'k:', 'linewidth', 2);hold off
    ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20)
end
legend({'APriD','MSA','CSA1','CSA2','PDSG-adp' },'Location','southwest','interpreter','latex')
xlabel('iteration','FontSize',20,'interpreter','latex')

subplot(2,3,2)
loglog(ks,out_scenario_ApriD.f1s_avgx_avg, 'b-', 'linewidth', 2);hold on
loglog(ks,out_scenario_MSA.f1s_avgx_avg, 'g--', 'linewidth', 2);hold on
loglog(ks,out_scenario_CSA.f1s_avgx_avg, 'r-.', 'linewidth', 2);hold on
loglog(ks,out_scenario_CSA.f1s_meanx_avg, 'r-', 'linewidth', 2);hold on
% loglog(ks,out_scenario_pdsg.f1s_avgx_avg, 'k--', 'linewidth', 2);hold on
loglog(ks,out_scenario_ada_pdsg.f1s_avgx_avg, 'k-', 'linewidth', 2);hold on
% loglog(ks,out_scenario_ada_pdsg1.f1s_avgx_avg, 'k-.', 'linewidth', 2);hold on
% loglog(ks,out_scenario_ada_pdsg2.f1s_avgx_avg, 'k:', 'linewidth', 2);hold on
% loglog(ks,out_scenario_ada_pdsg3.f1s_avgx_avg, 'c:', 'linewidth', 2);hold off
xlabel('iteration','FontSize',20,'interpreter','latex')
%ylabel('average $[f1(\mbox{x})]_+$','interpreter','latex','FontSize',20)
ylabel('average constraint violation','interpreter','latex','FontSize',20) 

subplot(2,3,3)
loglog(ks,out_scenario_ApriD.f1s_avgx_max, 'b-', 'linewidth', 2);hold on
loglog(ks,out_scenario_MSA.f1s_avgx_max, 'g--', 'linewidth', 2);hold on
loglog(ks,out_scenario_CSA.f1s_avgx_max, 'r-.', 'linewidth', 2);hold on
loglog(ks,out_scenario_CSA.f1s_meanx_max, 'r-', 'linewidth', 2);hold on
% loglog(ks,out_scenario_pdsg.f1s_avgx_max, 'k--', 'linewidth', 2);hold on
loglog(ks,out_scenario_ada_pdsg.f1s_avgx_max, 'k-', 'linewidth', 2);hold off
xlabel('iteration','FontSize',20,'interpreter','latex') 
ylabel('max constraint violation','interpreter','latex','FontSize',20)

subplot(2,3,4)
if n==10
    loglog(time_scenario_ApriD/ks(end)*ks,abs(out_scenario_ApriD.f0s_avgx-cvx_optval), 'b-', 'linewidth', 2);hold on
    loglog(time_scenario_MSA/ks(end)*ks,abs(out_scenario_MSA.f0s_avgx-cvx_optval), 'g--', 'linewidth', 2);hold on
    loglog(time_scenario_CSA/ks(end)*ks,abs(out_scenario_CSA.f0s_avgx-cvx_optval), 'r-.', 'linewidth', 2);hold on
    loglog(time_scenario_CSA/ks(end)*ks,abs(out_scenario_CSA.f0s_meanx-cvx_optval), 'r-', 'linewidth', 2);hold on 
    loglog(time_scenario_MSA/ks(end)*ks,abs(out_scenario_ada_pdsg.f0s_avgx-cvx_optval), 'k-', 'linewidth', 2);hold off
    ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20)
end
if n==200
    optval_estimate = min([
        min(out_scenario_ApriD.f0s_avgx(out_scenario_ApriD.f1s_avgx_max<=0)),...
        min(out_scenario_MSA.f0s_avgx(out_scenario_MSA.f1s_avgx_max<=0)),...
        min(out_scenario_CSA.f0s_avgx(out_scenario_CSA.f1s_avgx_max<=0)),... 
        min(out_scenario_ada_pdsg.f0s_avgx(out_scenario_ada_pdsg.f1s_avgx_max<=0))...
        ]);
    loglog(time_scenario_ApriD/ks(end)*ks,out_scenario_ApriD.f0s_avgx-optval_estimate, 'b-', 'linewidth', 2);hold on
    loglog(time_scenario_MSA/ks(end)*ks,out_scenario_MSA.f0s_avgx-optval_estimate, 'g--', 'linewidth', 2);hold on
    loglog(time_scenario_CSA/ks(end)*ks,out_scenario_CSA.f0s_avgx-optval_estimate, 'r-.', 'linewidth', 2);hold on
    loglog(time_scenario_CSA/ks(end)*ks,out_scenario_CSA.f0s_meanx-optval_estimate, 'r-', 'linewidth', 2);hold on 
    loglog(time_scenario_ada_pdsg/ks(end)*ks,out_scenario_ada_pdsg.f0s_avgx-optval_estimate, 'k-', 'linewidth', 2);hold off
    ylabel('$|f_0(\mathbf{x})-f_0(\mathbf{x}_{\mathbf{opt}})|$','interpreter','latex','FontSize',20)
end 
xlabel('time','FontSize',20,'interpreter','latex')
xlim([time_max/ks(end)*ks(1), time_max])

subplot(2,3,5)
loglog(time_scenario_ApriD/ks(end)*ks,out_scenario_ApriD.f1s_avgx_avg, 'b-', 'linewidth', 2);hold on
loglog(time_scenario_MSA/ks(end)*ks,out_scenario_MSA.f1s_avgx_avg, 'g--', 'linewidth', 2);hold on
loglog(time_scenario_CSA/ks(end)*ks,out_scenario_CSA.f1s_avgx_avg, 'r-.', 'linewidth', 2);hold on
loglog(time_scenario_CSA/ks(end)*ks,out_scenario_CSA.f1s_meanx_avg, 'r-', 'linewidth', 2);hold on 
loglog(time_scenario_ada_pdsg/ks(end)*ks,out_scenario_ada_pdsg.f1s_avgx_avg, 'k-', 'linewidth', 2);hold off 
xlabel('time','FontSize',20,'interpreter','latex')
%ylabel('average $[f1(\mbox{x})]_+$','interpreter','latex','FontSize',20)
ylabel('average constraint violation','interpreter','latex','FontSize',20) 
xlim([time_max/ks(end)*ks(1), time_max])

subplot(2,3,6)
loglog(time_scenario_ApriD/ks(end)*ks,out_scenario_ApriD.f1s_avgx_max, 'b-', 'linewidth', 2);hold on
loglog(time_scenario_MSA/ks(end)*ks,out_scenario_MSA.f1s_avgx_max, 'g--', 'linewidth', 2);hold on
loglog(time_scenario_CSA/ks(end)*ks,out_scenario_CSA.f1s_avgx_max, 'r-.', 'linewidth', 2);hold on
loglog(time_scenario_CSA/ks(end)*ks,out_scenario_CSA.f1s_meanx_max, 'r-', 'linewidth', 2);hold on 
loglog(time_scenario_ada_pdsg/ks(end)*ks,out_scenario_ada_pdsg.f1s_avgx_max, 'k-', 'linewidth', 2);hold off
xlabel('time','FontSize',20,'interpreter','latex')
%ylabel('maximum $[f1(\mbox{x})]_+$','interpreter','latex','FontSize',20)
ylabel('max constraint violation','interpreter','latex','FontSize',20)
xlim([time_max/ks(end)*ks(1), time_max])

sgtitle(['n = ' num2str(n) ', p = ' num2str(p)])
print(gcf,'-dpdf', savename);
print(gcf,'-depsc', savename);
