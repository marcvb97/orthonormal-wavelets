%% Analyze the results
%%
load KODAK_OW.mat
[av_OW, level_OW, kmax_OW] = analyze(RESULTS_GLOBAL, 3);

%%
load KODAK_OW2.mat
[av_OW2, level_OW2, kmax_OW2] = analyze(RESULTS_GLOBAL, 3);

%%
load KODAK_VP.mat
[av_VP, level_VP, kmax_VP] = analyze(RESULTS_GLOBAL, 3);

%%
load KODAK_db2.mat
[av_db2, level_db2, kmax_db2] = analyze(RESULTS_GLOBAL, 2);

%%
load KODAK_bior35.mat
[av_bior35, level_bior35, kmax_bior35] = analyze(RESULTS_GLOBAL, 2);

%%
load PEXELS300_OW_par.mat
[av_OW, level_OW, kmax_OW] = analyze(RESULTS_GLOBAL, 3);

%%
load PEXELS300_OW2_par.mat
[av_OW2, level_OW2, kmax_OW2] = analyze(RESULTS_GLOBAL, 3);

%%
load PEXELS300_VP_par.mat
[av_VP, level_VP, kmax_VP] = analyze(RESULTS_GLOBAL, 3);

%%
load PEXELS300_db2.mat
[av_db2, level_db2, kmax_db2] = analyze(RESULTS_GLOBAL, 2);

%%
load PEXELS300_bior35.mat
[av_bior35, level_bior35, kmax_bior35] = analyze(RESULTS_GLOBAL, 2);

%%
load NY96_OW_par.mat
[av_OW, level_OW, kmax_OW] = analyze(RESULTS_GLOBAL, 3);

%%
load NY96_OW2_par.mat
[av_OW2, level_OW2, kmax_OW2] = analyze(RESULTS_GLOBAL, 3);

%%
load NY96_VP_par.mat
[av_VP, level_VP, kmax_VP] = analyze(RESULTS_GLOBAL, 3);

%%
load NY96_db2.mat
[av_db2, level_db2, kmax_db2] = analyze(RESULTS_GLOBAL, 2);

%%
load NY96_bior35.mat
[av_bior35, level_bior35, kmax_bior35] = analyze(RESULTS_GLOBAL, 2);

%%
load 13US_OW.mat
[av_OW, level_OW, kmax_OW] = analyze(RESULTS_GLOBAL, 3);

%%
load 13US_OW2.mat
[av_OW2, level_OW2, kmax_OW2] = analyze(RESULTS_GLOBAL, 3);

%%
load 13US_VP.mat
[av_VP, level_VP, kmax_VP] = analyze(RESULTS_GLOBAL, 3);

%%
load 13US_db2.mat
[av_db2, level_db2, kmax_db2] = analyze(RESULTS_GLOBAL, 2);

%%
load 13US_bior35.mat
[av_bior35, level_bior35, kmax_bior35] = analyze(RESULTS_GLOBAL, 2);

%%
load URBAN100_OW.mat
[av_OW, level_OW, kmax_OW] = analyze(RESULTS_GLOBAL, 3);

%%
load URBAN100_OW2.mat
[av_OW2, level_OW2, kmax_OW2] = analyze(RESULTS_GLOBAL, 3);

%%
load URBAN100_VP.mat
[av_VP, level_VP, kmax_VP] = analyze(RESULTS_GLOBAL, 3);

%%
load URBAN100_db2.mat
[av_db2, level_db2, kmax_db2] = analyze(RESULTS_GLOBAL, 2);

%%
load URBAN100_bior35.mat
[av_bior35, level_bior35, kmax_bior35] = analyze(RESULTS_GLOBAL, 2);

%%
load NY17_OW_par.mat
[av_OW, level_OW, kmax_OW] = analyze(RESULTS_GLOBAL, 3);

%%
load NY17_OW2_par.mat
[av_OW2, level_OW2, kmax_OW2] = analyze(RESULTS_GLOBAL, 3);

%%
load NY17_VP_par.mat
[av_VP, level_VP, kmax_VP] = analyze(RESULTS_GLOBAL, 3);

%%
load NY17_db2.mat
[av_db2, level_db2, kmax_db2] = analyze(RESULTS_GLOBAL, 2);

%%
load NY17_bior35.mat
[av_bior35, level_bior35, kmax_bior35] = analyze(RESULTS_GLOBAL, 2);

%%
clf
set(0, 'DefaultAxesFontSize', 16);    % Axes labels & ticks
set(0, 'DefaultTextFontSize', 16);    % Text & titles
set(0, 'DefaultLineLineWidth', 2);  % All future plots will use linewidth = 2
set(0, 'DefaultLineMarkerSize', 10);
set(0, 'DefaultLegendFontSize', 16);

%%
clf
semilogx(av_OW(:,1),av_OW(:,2),'-o')
hold on
semilogx(av_OW2(:,1),av_OW2(:,2),'-^')
semilogx(av_VP(:,1),av_VP(:,2),'-x')
semilogx(av_db2(:,1),av_db2(:,2),'-+')
semilogx(av_bior35(:,1),av_bior35(:,2),'-v')
title('KODAK','Interpreter','latex')
xlabel("percentage",'Interpreter', 'latex')
ylabel("average PSNR",'Interpreter', 'latex')
legend("OVP","OVP2","IVP","DB2","BIOR3.5","Location","northwest",'Interpreter', 'latex')
saveas(1,'figs/fig_KODAK_PSNR.eps','epsc')

%%
clf
semilogx(av_OW(:,1),av_OW(:,3),'-o')
hold on
semilogx(av_OW2(:,1),av_OW2(:,3),'-^')
semilogx(av_VP(:,1),av_VP(:,3),'-x')
semilogx(av_db2(:,1),av_db2(:,3),'-+')
semilogx(av_bior35(:,1),av_bior35(:,3),'-v')
title('KODAK','Interpreter','latex')
xlabel("percentage",'Interpreter', 'latex')
ylabel("average SSSIM",'Interpreter', 'latex')
legend("OVP","OVP2","IVP","DB2","BIOR3.5","Location","northwest",'Interpreter', 'latex')
saveas(1,'figs/fig_KODAK_SSIM.eps','epsc')

%%
clf
semilogx(av_OW(:,1),av_OW(:,6),'-o')
hold on
semilogx(av_OW2(:,1),av_OW2(:,6),'-^')
semilogx(av_VP(:,1),av_VP(:,6),'-x')
semilogx(av_db2(:,1),av_db2(:,6),'-+')
semilogx(av_bior35(:,1),av_bior35(:,6),'-v')
title('KODAK','Interpreter','latex')
xlabel("percentage",'Interpreter', 'latex')
ylabel("average level",'Interpreter', 'latex')
legend("OVP","OVP2","IVP","DB2","BIOR3.5","Location","northwest",'Interpreter', 'latex')
saveas(1,'figs/fig_KODAK_level.eps','epsc')

%%
clf
semilogx(av_OW(:,1),av_OW(:,5),'-o')
hold on
semilogx(av_OW2(:,1),av_OW2(:,5),'-^')
semilogx(av_VP(:,1),av_VP(:,5),'-x')
title('KODAK','Interpreter','latex')
xlabel("percentage",'Interpreter', 'latex')
ylabel("average $\theta$",'Interpreter', 'latex')
legend("OVP","OVP2","IVP","Location","northwest",'Interpreter', 'latex')
saveas(1,'figs/fig_KODAK_theta.eps','epsc')

%%
