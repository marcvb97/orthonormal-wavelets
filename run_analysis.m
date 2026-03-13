%% Analyze the results

function run_analysis(dataset_name)

%%
get_filename(dataset_name, 'OW')
load(get_filename(dataset_name, 'OW'), 'RESULTS_GLOBAL');
[av_OW, level_OW, kmax_OW] = analyze(RESULTS_GLOBAL, 3);

%%
load(get_filename(dataset_name, 'OW2'), 'RESULTS_GLOBAL');
[av_OW2, level_OW2, kmax_OW2] = analyze(RESULTS_GLOBAL, 3);

%%
load(get_filename(dataset_name, 'VP'), 'RESULTS_GLOBAL');
[av_VP, level_VP, kmax_VP] = analyze(RESULTS_GLOBAL, 3);

%%
load(get_filename(dataset_name, 'db2'), 'RESULTS_GLOBAL');
[av_db2, level_db2, kmax_db2] = analyze(RESULTS_GLOBAL, 2);

%%
load(get_filename(dataset_name, 'bior35'), 'RESULTS_GLOBAL');
[av_bior35, level_bior35, kmax_bior35] = analyze(RESULTS_GLOBAL, 2);

%%
clf
set(0, 'DefaultAxesFontSize', 16);
set(0, 'DefaultTextFontSize', 16);
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultLineMarkerSize', 10);
set(0, 'DefaultLegendFontSize', 16);

%%
clf
semilogx(av_OW(:,1)*100,av_OW(:,2),'-o')
hold on
semilogx(av_OW2(:,1)*100,av_OW2(:,2),'-^')
semilogx(av_VP(:,1)*100,av_VP(:,2),'-x')
semilogx(av_db2(:,1)*100,av_db2(:,2),'-+')
semilogx(av_bior35(:,1)*100,av_bior35(:,2),'-v')
title(dataset_name,'Interpreter','latex')
xlabel("percentage",'Interpreter', 'latex')
ylabel("average PSNR",'Interpreter', 'latex')
legend("OVP","OVP2","IVP","DB2","BIOR3.5","Location","northwest",'Interpreter', 'latex')
saveas(1, sprintf('figs/fig_%s_PSNR.eps', dataset_name), 'epsc')

%%
clf
semilogx(av_OW(:,1)*100,av_OW(:,3),'-o')
hold on
semilogx(av_OW2(:,1)*100,av_OW2(:,3),'-^')
semilogx(av_VP(:,1)*100,av_VP(:,3),'-x')
semilogx(av_db2(:,1)*100,av_db2(:,3),'-+')
semilogx(av_bior35(:,1)*100,av_bior35(:,3),'-v')
title(dataset_name,'Interpreter','latex')
xlabel("percentage",'Interpreter', 'latex')
ylabel("average SSSIM",'Interpreter', 'latex')
legend("OVP","OVP2","IVP","DB2","BIOR3.5","Location","northwest",'Interpreter', 'latex')
saveas(1, sprintf('figs/fig_%s_SSIM.eps', dataset_name), 'epsc')

%%
clf
semilogx(av_OW(:,1)*100,av_OW(:,6),'-o')
hold on
semilogx(av_OW2(:,1)*100,av_OW2(:,6),'-^')
semilogx(av_VP(:,1)*100,av_VP(:,6),'-x')
semilogx(av_db2(:,1)*100,av_db2(:,6),'-+')
semilogx(av_bior35(:,1)*100,av_bior35(:,6),'-v')
title(dataset_name,'Interpreter','latex')
xlabel("percentage",'Interpreter', 'latex')
ylabel("average level",'Interpreter', 'latex')
legend("OVP","OVP2","IVP","DB2","BIOR3.5","Location","northwest",'Interpreter', 'latex')
saveas(1, sprintf('figs/fig_%s_level.eps', dataset_name), 'epsc')

%%
clf
semilogx(av_OW(:,1)*100,av_OW(:,5),'-o')
hold on
semilogx(av_OW2(:,1)*100,av_OW2(:,5),'-^')
semilogx(av_VP(:,1)*100,av_VP(:,5),'-x')
title(dataset_name,'Interpreter','latex')
xlabel("percentage",'Interpreter', 'latex')
ylabel("average $\theta$",'Interpreter', 'latex')
legend("OVP","OVP2","IVP","Location","northwest",'Interpreter', 'latex')
saveas(1, sprintf('figs/fig_%s_theta.eps', dataset_name), 'epsc')

end

% -------------------------------------------------------------------------
function filename = get_filename(dataset_name, suffix)
    par_datasets = {'PEXELS300', 'NY96', 'NY17'};
    par_suffixes = {'OW', 'OW2', 'VP'};

    if ismember(dataset_name, par_datasets) && ismember(suffix, par_suffixes)
        filename = sprintf('mat/%s_%s_par.mat', dataset_name, suffix);
    else
        filename = sprintf('mat/%s_%s.mat', dataset_name, suffix);
    end
end