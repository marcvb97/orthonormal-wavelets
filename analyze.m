function [av, level, kmax] = analyze(res, factor)
RESULTS_best = res{1}.RESULTS_best;
[m,~] =size(RESULTS_best);
av = RESULTS_best;
max_size = max(res{1}.size);
kmax = log(max_size)/log(factor) * ones(m,1);
level = res{1}.RESULTS_best(:,end);
for i = 2:length(res)
    RESULTS_best = res{i}.RESULTS_best;
    [m,~] =size(RESULTS_best);
    av = av + RESULTS_best;
    max_size = max(res{i}.size);
    kmax = [kmax; log(max_size)/log(factor) * ones(m,1)];
    level = [level; RESULTS_best(:,end)];
end
av = av / length(res);