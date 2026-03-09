%% timings for the 1D transformations

powmax = 12;
tt = zeros(powmax,4);
theta = 0.5;
fsca = sqrt(3); fwav = sqrt(3/2);
kmax = 1000;
for pow = 5:powmax
    pow
    n = 2^pow;
    v = rand(n,1);
    tic
    for k = 1:kmax
        [anew, bnew, ~] = FT1step_m(v, theta);
    end
    tt(pow,1) = tt(pow,1) + toc;
    tic
    for k = 1:kmax
        [anew, bnew, nplus] = FWT1step_m_VP(v, theta,fsca,fwav);
    end
    tt(pow,2) = tt(pow,2) + toc;
    tic
    for k = 1:kmax
        [IdecR,SdecR]=wavedec(v,2,"db2");
    end
    tt(pow,3) = tt(pow,3) + toc;
    tic
    for k = 1:kmax
        [IdecR,SdecR]=wavedec(v,2,"bior3.5");
    end
    tt(pow,4) = tt(pow,4) + toc;
end

tt = tt / kmax;

%%
semilogy(tt,'-o')
xlabel("power of 2")
ylabel("time (seconds)")
legend("OW","VP","DB2","BIOR3.5")
title("timings 1D")

%%
semilogy(tt(:,1)./tt(:,3),'o-')
hold on
semilogy(tt(:,2)./tt(:,3),'+-')
semilogy(tt(:,4)./tt(:,3),'x-')
xlabel("power of 2")
ylabel("factor")
legend("OW","VP","BIOR3.5")
title("timings 1D relative to DB2")
