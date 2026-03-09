function [TT, kmax] = analyze_timings(RESULTS_GLOBAL,factor, kmin)

for pow = 5:9
    kmax(pow) = floor(log(2^pow)/log(factor) + 1.0e-12);
    for k = kmin:kmin   % kmax(pow)
        T = [];
        for j = 1:length(RESULTS_GLOBAL)
            times = RESULTS_GLOBAL{j}.TIMES;
            for l = 1:size(times,1)
                if (RESULTS_GLOBAL{j}.size(1) == 2^pow) && (times(l,1) == k)
                    T = [T, times(l,4)+times(l,5)+times(l,6)];
                end
            end
        end
        TT(pow,k) = sum(T) / length(T);
    end
end

