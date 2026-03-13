function [RESULTS_GLOBAL] = imagecompression_color_other_PSNR_SSIM(filepath)

% AIM: Image compression using custom wavelets.
% The best result is found by varying theta (see vtheta).
% Depending on the value of flag, all decomposition coefficients 
% (wavelet and scaling) are compressed or only the wavelet coefficients.
% Results are displayed for all possible decomposition steps (see Kmax).
% initial_transformation indicates if TR2D and ITR2D are used.
%
% MODIFIED: search now independently minimises MSE and maximises SSIM.

%% Settings
flag = 1;   % 0 = compress all coefficients, 1 = compress details only
threshold_type = 1; % 1 = separately, 2 = combined over the channels

vkeep  = [0.5, 0.25, 0.10, 0.05, 0.02, 0.01, 0.005];
% vkeep  = [0.5];
vtheta = [0.0];

wname='db2';

fsca=sqrt(3); fwav=sqrt(3/2);

%% Load and prepare image
testImage = filepath;
II = imread(testImage);
imshow(II)
R = double(II(:,:,1));
G = double(II(:,:,2));
B = double(II(:,:,3));

[Nor, Mor] = size(R);

%% Decomposition parameters
Kmax = floor(log(min(Nor,Mor))/log(2) + 1.0e-12);
Kmax = 6;
dmin = 3;

%% Main loop
% RESULTS columns: [keep, PSNR_mse, SSIM_mse, %retained_mse, theta_mse, kmax_mse,
%                         PSNR_ssim, SSIM_ssim, theta_ssim, kmax_ssim]
RESULTS = zeros((Kmax-4+1) * length(vkeep), 10);
TIMES   = zeros((Kmax-4+1) * length(vkeep) * length(vtheta), 7);
cont = 0;
cont_time = 0;

for kstep = 4:Kmax
    for keep = vkeep
        cont = cont + 1;
        cont_time = cont_time + 1;

        %% First theta: initialise both trackers
        theta = vtheta(1);

        tic;
        [mseVAL, ssimVAL, IoutR, IoutG, IoutB, zero_el, IcompR] = ...
            run_theta(R, G, B, kstep, keep, flag, threshold_type, Nor, Mor, wname);
        t_total = toc;

        % MSE tracker
        best_mse      = mseVAL;
        Rfin_mse      = IoutR;   Gfin_mse = IoutG;   Bfin_mse = IoutB;
        theta_fin_mse = theta;   kmax_mse = kstep;
        zero_el_mse   = zero_el;
        IcompR_mse    = IcompR;

        % SSIM tracker
        best_ssim      = ssimVAL;
        Rfin_ssim      = IoutR;   Gfin_ssim = IoutG;   Bfin_ssim = IoutB;
        theta_fin_ssim = theta;   kmax_ssim = kstep;
        zero_el_ssim   = zero_el;

        TIMES(cont_time, :) = [kstep, keep, theta, t_total, 0, 0, mseVAL];

        %% Search remaining theta values (none in this code by default,
        %  but kept symmetric with imagecompression_color_OW for extensibility)
        for theta = setdiff(vtheta, vtheta(1))
            cont_time = cont_time + 1;
            tic;
            [mseVAL, ssimVAL, IoutR, IoutG, IoutB, zero_el, IcompR] = ...
                run_theta(R, G, B, kstep, keep, flag, threshold_type, Nor, Mor, wname);
            t_total = toc;

            % Update MSE-optimal tracker
            if mseVAL < best_mse
                best_mse      = mseVAL;
                Rfin_mse      = IoutR;   Gfin_mse = IoutG;   Bfin_mse = IoutB;
                theta_fin_mse = theta;   kmax_mse = kstep;
                zero_el_mse   = zero_el;
                IcompR_mse    = IcompR;
            end

            % Update SSIM-optimal tracker
            if ssimVAL > best_ssim
                best_ssim      = ssimVAL;
                Rfin_ssim      = IoutR;   Gfin_ssim = IoutG;   Bfin_ssim = IoutB;
                theta_fin_ssim = theta;   kmax_ssim = kstep;
                zero_el_ssim   = zero_el;
            end

            TIMES(cont_time, :) = [kstep, keep, theta, t_total, 0, 0, mseVAL];
        end

        %% Quality metrics for the MSE-optimal result
        disp(filepath)
        psnr_mse = 10 * log10(255^2 / best_mse);
        ssim_mse = (ssim(R, Rfin_mse) + ssim(G, Gfin_mse) + ssim(B, Bfin_mse)) / 3;

        %% Quality metrics for the SSIM-optimal result
        DR = abs(double(R) - double(Rfin_ssim)).^2;
        DG = abs(double(G) - double(Gfin_ssim)).^2;
        DB = abs(double(B) - double(Bfin_ssim)).^2;
        mse_ssim  = sum(DR(:) + DG(:) + DB(:)) / (numel(R) + numel(G) + numel(B));
        psnr_ssim = 10 * log10(255^2 / mse_ssim);
        ssim_ssim = best_ssim;

        num = 3 * length(IcompR_mse);

        pct_retained_mse = (num - zero_el_mse) * 100 / (3*Nor*Mor);

        if pct_retained_mse < 70
            % Columns: [keep | PSNR_mse  SSIM_mse  %retained_mse  theta_mse  kmax_mse |
            %                  PSNR_ssim SSIM_ssim  theta_ssim     kmax_ssim]
            RESULTS(cont, :) = [keep, ...
                psnr_mse,  ssim_mse,  pct_retained_mse, theta_fin_mse, kmax_mse, ...
                psnr_ssim, ssim_ssim, theta_fin_ssim,   kmax_ssim]

            % Show both optimal reconstructions
            figure(2*cont - 1); imshow(uint8(cat(3, Rfin_mse, Gfin_mse, Bfin_mse)));
            title(sprintf('MSE-optimal  PSNR=%.2f dB', psnr_mse));
            axis off equal

            figure(2*cont); imshow(uint8(cat(3, Rfin_ssim, Gfin_ssim, Bfin_ssim)));
            title(sprintf('SSIM-optimal  SSIM=%.4f', ssim_ssim));
            axis off equal
        else
            disp("percentage retained >= 70%")
            cont      = cont - 1;
            cont_time = cont_time - 1;
        end

    end
end

%% Display results table
varNames = {'CR', ...
    'PSNR_{MSE}', 'SSIM_{MSE}', '%retained_{MSE}', 'theta_{MSE}', 'Decomp_{MSE}', ...
    'PSNR_{SSIM}', 'SSIM_{SSIM}', 'theta_{SSIM}', 'Decomp_{SSIM}'};

% Best by PSNR (MSE-optimal): maximise column 2 = PSNR_mse
RESULTS_best_psnr = zeros(length(vkeep), size(RESULTS, 2));
for i = 1:length(vkeep)
    I = find(abs(RESULTS(:,1) - vkeep(i)) < 1e-12 & RESULTS(:,2) ~= 0);
    [~, J] = max(RESULTS(I, 2));
    RESULTS_best_psnr(i, :) = RESULTS(I(J(1)), :);
end

% Best by SSIM (SSIM-optimal): maximise column 8 = SSIM_ssim
RESULTS_best_ssim = zeros(length(vkeep), size(RESULTS, 2));
for i = 1:length(vkeep)
    I = find(abs(RESULTS(:,1) - vkeep(i)) < 1e-12 & RESULTS(:,2) ~= 0);
    [~, J] = max(RESULTS(I, 8));
    RESULTS_best_ssim(i, :) = RESULTS(I(J(1)), :);
end

disp('***** BEST BY PSNR *****')
T_psnr = array2table(RESULTS_best_psnr, 'VariableNames', varNames);
disp(T_psnr)

disp('***** BEST BY SSIM *****')
T_ssim = array2table(RESULTS_best_ssim, 'VariableNames', varNames);
disp(T_ssim)

%% RESULTS_GLOBAL
RESULTS_GLOBAL = struct;
RESULTS_GLOBAL.name              = filepath;
RESULTS_GLOBAL.RESULTS           = RESULTS;
RESULTS_GLOBAL.RESULTS_best_psnr = RESULTS_best_psnr;
RESULTS_GLOBAL.RESULTS_best_ssim = RESULTS_best_ssim;
RESULTS_GLOBAL.TIMES             = TIMES;
RESULTS_GLOBAL.size              = [Nor, Mor];

end % end of main function


%% --- Local helper function ---
function [mseVAL, ssimVAL, IoutR_, IoutG_, IoutB_, zero_el_, IcompR_] = ...
        run_theta(R, G, B, kstep, keep, flag, threshold_type, Nor, Mor, wname)

    [IdecR, SdecR] = wavedec2(R, kstep, wname);
    [IdecG, SdecG] = wavedec2(G, kstep, wname);
    [IdecB, SdecB] = wavedec2(B, kstep, wname);

    num_keep = keep * Nor * Mor;
    keep2    = num_keep / length(IdecR);

    switch threshold_type
        case 1
            [IcompR_, zero_elR] = THR_other(IdecR, SdecR, keep2, flag);
            [IcompG,  zero_elG] = THR_other(IdecG, SdecG, keep2, flag);
            [IcompB,  zero_elB] = THR_other(IdecB, SdecB, keep2, flag);
            zero_el_ = zero_elR + zero_elG + zero_elB;
        case 2
            [~, zero_el_, IcompR_, IcompG, IcompB] = THR_color( ...
                IdecR, IdecG, IdecB, [], [], [], [], [], [], keep, [], [], flag);
    end

    IrecR = waverec2(IcompR_, SdecR, wname);
    IrecG = waverec2(IcompG,  SdecG, wname);
    IrecB = waverec2(IcompB,  SdecB, wname);

    IoutR_ = IrecR(1:Nor, 1:Mor);
    IoutG_ = IrecG(1:Nor, 1:Mor);
    IoutB_ = IrecB(1:Nor, 1:Mor);

    DR = abs(double(R) - double(IoutR_)).^2;
    DG = abs(double(G) - double(IoutG_)).^2;
    DB = abs(double(B) - double(IoutB_)).^2;
    mseVAL  = sum(DR(:) + DG(:) + DB(:)) / (numel(R) + numel(G) + numel(B));
    ssimVAL = (ssim(R, IoutR_) + ssim(G, IoutG_) + ssim(B, IoutB_)) / 3;

end