%% PSNR in function of theta and level
T = RESULTS_GLOBAL{1}.TIMES;
T = reshape(T(:,end),9,3*7);
T = 10 * log10(255^2 ./ T);
ttheta = 0.1:0.1:0.9; ttheta = ttheta.';
plot(ttheta, T,'-o')
xlabel("$\theta$",'Interpreter','latex')
ylabel("PSNR",'Interpreter','latex')
legend('50\%','25\%','10\%','5\%','2\%','1\%','0.5\%','Interpreter','latex')
saveas(1,'fig_PSNR_theta_level.eps','epsc')
