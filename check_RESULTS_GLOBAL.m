%% what levels and theta values are present in RESULTS_GLOBAL
level_min = Inf;
level_max = 0;
ttheta = [];
for i = 1:length(RESULTS_GLOBAL)
    RES = RESULTS_GLOBAL{i}.RESULTS;
    level_min = min([level_min;RES(:,6)]);
    level_max = max([level_max;RES(:,6)]);
    ttheta = [ttheta; RES(:,5)];
end
length(RESULTS_GLOBAL)
[level_min, level_max]
hist(ttheta)