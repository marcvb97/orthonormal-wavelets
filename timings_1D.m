%% timings for the 1D transformations

powmax = 12;
tt = zeros(powmax,5);
theta = 0.5;
fsca = sqrt(3); fwav = sqrt(3/2);
kmax = 1000;
FT1step_m(rand(3*4,1), 0.5);          % dummy run
FWT1step_m_VP(rand(3*4,1), 0.5, 1, 1);
dct(rand(10,1))
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
        % vt = dct(v);
    end
    tt(pow,3) = tt(pow,3) + toc;
    tic
    for k = 1:kmax
        [IdecR,SdecR]=wavedec(v,2,"bior3.5");
    end
    tt(pow,4) = tt(pow,4) + toc;
    v3 = rand(3*n,kmax);
    % tic
    % for k = 1:kmax
    %     v = dct(v3(:,k));
    % end
    % tt(pow,5) = tt(pow,5) + toc;
end

tt = tt / kmax;

%%
clf
set(0, 'DefaultAxesFontSize', 16);    % Axes labels & ticks
set(0, 'DefaultTextFontSize', 16);    % Text & titles
set(0, 'DefaultLineLineWidth', 2);  % All future plots will use linewidth = 2
set(0, 'DefaultLineMarkerSize', 10);
set(0, 'DefaultLegendFontSize', 16);

%%
clf
semilogy(5:12,tt(5:end,1),'-o')
hold on
semilogy(5:12,tt(5:end,2),'-+')
semilogy(5:12,tt(5:end,3),'-x')
semilogy(5:12,tt(5:end,4),'-v')
xlabel("power of 2",'Interpreter','latex')
ylabel("time (seconds)",'Interpreter','latex')
legend("OVP","IVP","DB2","BIOR3.5",'Interpreter','latex','Location','NW')
title("timings 1D",'Interpreter','latex')
saveas(1,['figs/fig_timings1D_absolute'],'epsc')

%%
clf
plot(tt(:,1)./tt(:,3),'o-')
hold on
plot(tt(:,2)./tt(:,3),'+-')
plot(tt(:,4)./tt(:,3),'x-')
xlabel("power of 2",'Interpreter','latex')
ylabel("factor",'Interpreter','latex')
legend("OVP","IVP","BIOR3.5",'Interpreter','latex','Location','NW')
title("timings 1D relative to DB2",'Interpreter','latex')
saveas(1,'figs/fig_timings1D_relative','epsc')

