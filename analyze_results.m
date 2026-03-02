%% Analyze the results
%%
load KODAK_OW.mat
[av_OW, level_OW, kmax_OW] = analyze(RESULTS_GLOBAL, 3);

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
load NY96_OW.mat
[av_OW, level_OW, kmax_OW] = analyze(RESULTS_GLOBAL, 3);

%%
load NY96_VP.mat
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
load URBAN100_VP.mat
[av_VP, level_VP, kmax_VP] = analyze(RESULTS_GLOBAL, 3);

%%
load URBAN100_db2.mat
[av_db2, level_db2, kmax_db2] = analyze(RESULTS_GLOBAL, 2);

%%
load URBAN100_bior35.mat
[av_bior35, level_bior35, kmax_bior35] = analyze(RESULTS_GLOBAL, 2);


%%
semilogx(av_OW(:,1),av_OW(:,2),'-o')
hold on
% semilogx(av_VP(:,1),av_VP(:,2),'-x')
semilogx(av_db2(:,1),av_db2(:,2),'-+')
semilogx(av_bior35(:,1),av_bior35(:,2),'-v')

