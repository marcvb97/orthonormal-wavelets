%% test FT1step and IFT1step
theta = 0.5;
a = rand(1000,1);
[anew, bnew, nplus] = FT1step_m(a, theta);
[aa, nminus] = IFT1step_m(anew, bnew, theta);
norm(aa(1:end-nplus)-a)