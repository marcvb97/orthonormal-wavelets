%% analyze timings
load RANDOM_OW.mat

kmin = 4; factor = 3;
[TT, kmax] = analyze_timings(RESULTS_GLOBAL,factor, kmin);

semilogy(TT,'-o')
hold on

%%
load RANDOM_OW_par.mat

kmin = 4; factor = 3;
[TT, kmax] = analyze_timings(RESULTS_GLOBAL,factor, kmin);

semilogy(TT,'-x')

%%
load RANDOM_OW_par2.mat

kmin = 4; factor = 3;
[TT, kmax] = analyze_timings(RESULTS_GLOBAL,factor, kmin);

semilogy(TT,'-+')

%%
load RANDOM_db2.mat

kmin = 4; factor = 2;
[TT, kmax] = analyze_timings(RESULTS_GLOBAL,factor, kmin);

semilogy(TT,'-^')

%%
load RANDOM_bior35.mat

kmin = 4; factor = 2;
[TT, kmax] = analyze_timings(RESULTS_GLOBAL,factor, kmin);

semilogy(TT,'-v')

xlabel("power of 2")
ylabel("execution time (seconds)")
% legend("OW","OW2","VP","DB2","BIOR3.5")
