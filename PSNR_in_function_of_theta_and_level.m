%% PSNR in function of theta and level
T = RESULTS_GLOBAL{1}.TIMES;
T = reshape(T,9,3*7);
T = 10 * log10(255^2 ./ T);
ttheta = 0.1:0.1:0.9; ttheta = ttheta.';
plot(ttheta, T,'-o')
xlabel("\theta")
ylabel("PSNR")