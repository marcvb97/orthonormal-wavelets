function [RESULTS_GLOBAL] = imagecompression_color_OW_PSNR_SSIM(filepath)

% AIM: Image compression using custom wavelets.
% The best result is found by varying theta (see vtheta).
% Depending on the value of flag, all decomposition coefficients 
% (wavelet and scaling) are compressed or only the wavelet coefficients.
% Results are displayed for all possible decomposition steps (see Kmax).
% initial_transformation indicates if TR2D and ITR2D are used.
%
% MODIFIED: theta search now independently minimises MSE and maximises SSIM.

%% Settings
flag = 1;   % 0 = compress all coefficients, 1 = compress details only
threshold_type = 1; % 1 = separately, 2 = combined over the channels
initial_transformation = true;

vkeep  = [0.5, 0.25, 0.10, 0.05, 0.02, 0.01, 0.005];
% vkeep  = [0.5];
vtheta = [0.1:0.1:0.9];
% vtheta = [0.5];

%% Load and prepare image
testImage = filepath;
II = imread(testImage);
imshow(II)
R = double(II(:,:,1));
G = double(II(:,:,2));
B = double(II(:,:,3));

[Nor, Mor] = size(R);

%% Decomposition parameters
Kmax = floor(log(min(Nor, Mor))/log(3) + 1.0e-12);
Kmax = 5;
dmin = 3;

%% Main loop
% RESULTS columns: [keep, PSNR_mse, SSIM_mse, %retained_mse, theta_mse, kmax_mse,
%                         PSNR_ssim, SSIM_ssim, theta_ssim, kmax_ssim]
RESULTS = zeros((Kmax-3+1) * length(vkeep), 10);
TIMES   = zeros((Kmax-3+1) * length(vkeep) * length(vtheta), 7);
cont = 0;
cont_time = 0;

for kstep = 3:Kmax
    for keep = vkeep
        cont = cont + 1;
        cont_time = cont_time + 1;

        %% First theta: initialise both trackers
        theta = vtheta(1);
        tic;
        [mseVAL, ssimVAL, IoutR, IoutG, IoutB, lrow, lcol, zero_el] = ...
            run_theta(R, G, B, theta, kstep, dmin, keep, flag, threshold_type, ...
                      initial_transformation, Nor, Mor);
        t_total = toc;

        % MSE tracker
        best_mse      = mseVAL;
        Rfin_mse      = IoutR;   Gfin_mse = IoutG;   Bfin_mse = IoutB;
        theta_fin_mse = theta;   kmax_mse = max(lrow, lcol);
        zero_el_mse   = zero_el;

        % SSIM tracker
        best_ssim      = ssimVAL;
        Rfin_ssim      = IoutR;   Gfin_ssim = IoutG;   Bfin_ssim = IoutB;
        theta_fin_ssim = theta;   kmax_ssim = max(lrow, lcol);
        zero_el_ssim   = zero_el;

        TIMES(cont_time, :) = [kstep, keep, theta, t_total, 0, 0, mseVAL];

        %% Search remaining theta values
        for theta = setdiff(vtheta, vtheta(1))
            cont_time = cont_time + 1;
            tic;
            [mseVAL, ssimVAL, IoutR, IoutG, IoutB, lrow, lcol, zero_el] = ...
                run_theta(R, G, B, theta, kstep, dmin, keep, flag, threshold_type, ...
                          initial_transformation, Nor, Mor);
            t_total = toc;

            % Update MSE-optimal tracker
            if mseVAL < best_mse
                best_mse      = mseVAL;
                Rfin_mse      = IoutR;   Gfin_mse = IoutG;   Bfin_mse = IoutB;
                theta_fin_mse = theta;   kmax_mse = max(lrow, lcol);
                zero_el_mse   = zero_el;
            end

            % Update SSIM-optimal tracker
            if ssimVAL > best_ssim
                best_ssim      = ssimVAL;
                Rfin_ssim      = IoutR;   Gfin_ssim = IoutG;   Bfin_ssim = IoutB;
                theta_fin_ssim = theta;   kmax_ssim = max(lrow, lcol);
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

        % Retrieve IscalR/IwavR for num calculation via one final DEC call
        if initial_transformation
            AR = TR2D(R, theta_fin_mse);
        else
            AR = R;
        end
        [~, IscalR, IwavR, ~, ~] = DEC(AR, theta_fin_mse, kstep, dmin);
        num = 3 * length([IscalR; IwavR]);

        % Columns: [keep | PSNR_mse  SSIM_mse  %retained_mse  theta_mse  kmax_mse |
        %                  PSNR_ssim SSIM_ssim  theta_ssim     kmax_ssim]
        RESULTS(cont, :) = [keep, ...
            psnr_mse,  ssim_mse,  (num - zero_el_mse) * 100 / (3*Nor*Mor), theta_fin_mse,  kmax_mse, ...
            psnr_ssim, ssim_ssim, theta_fin_ssim,                           kmax_ssim]

        % Show both optimal reconstructions
        figure(2*cont - 1); imshow(uint8(cat(3, Rfin_mse,  Gfin_mse,  Bfin_mse)));
        title(sprintf('MSE-optimal  \\theta=%.1f  PSNR=%.2f dB', theta_fin_mse, psnr_mse));
        axis off equal

        figure(2*cont);     imshow(uint8(cat(3, Rfin_ssim, Gfin_ssim, Bfin_ssim)));
        title(sprintf('SSIM-optimal \\theta=%.1f  SSIM=%.4f', theta_fin_ssim, ssim_ssim));
        axis off equal
    end
end

%% Display results table
varNames = {'CR', ...
    'PSNR_{MSE}', 'SSIM_{MSE}', '%retained_{MSE}', 'theta_{MSE}', 'Decomp_{MSE}', ...
    'PSNR_{SSIM}', 'SSIM_{SSIM}', 'theta_{SSIM}', 'Decomp_{SSIM}'};

% Best by PSNR (MSE-optimal theta): maximise column 2 = PSNR_mse
RESULTS_best_psnr = zeros(length(vkeep), size(RESULTS, 2));
for i = 1:length(vkeep)
    I = find(abs(RESULTS(:,1) - vkeep(i)) < 1e-12 & RESULTS(:,2) ~= 0);
    [~, J] = max(RESULTS(I, 2));
    RESULTS_best_psnr(i, :) = RESULTS(I(J(1)), :);
end

% Best by SSIM (SSIM-optimal theta): maximise column 8 = SSIM_ssim
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
function [mseVAL, ssimVAL, IoutR_, IoutG_, IoutB_, lrow_, lcol_, zero_el_] = ...
        run_theta(R, G, B, theta, kstep, dmin, keep, flag, threshold_type, ...
                  initial_transformation, Nor, Mor)

    if initial_transformation
        AR = TR2D(R, theta);  AG = TR2D(G, theta);  AB = TR2D(B, theta);
    else
        AR = R;  AG = G;  AB = B;
    end

    [IdecR, IscalR, IwavR, lrow_, lcol_] = DEC(AR, theta, kstep, dmin);
    [IdecG, IscalG, IwavG,      ~,      ~] = DEC(AG, theta, kstep, dmin);
    [IdecB, IscalB, IwavB,      ~,      ~] = DEC(AB, theta, kstep, dmin);

    switch threshold_type
        case 1
            [~, zero_elR, IcompR] = THR(IdecR, IscalR, IwavR, keep, lrow_, lcol_, flag);
            [~, zero_elG, IcompG] = THR(IdecG, IscalG, IwavG, keep, lrow_, lcol_, flag);
            [~, zero_elB, IcompB] = THR(IdecB, IscalB, IwavB, keep, lrow_, lcol_, flag);
            zero_el_ = zero_elR + zero_elG + zero_elB;
        case 2
            [~, zero_el_, IcompR, IcompG, IcompB] = THR_color( ...
                IdecR, IdecG, IdecB, IscalR, IscalG, IscalB, ...
                IwavR, IwavG, IwavB, keep, lrow_, lcol_, flag);
    end

    IrecR = IFWTmatrix_m(IcompR, lrow_, lcol_, theta);
    IrecG = IFWTmatrix_m(IcompG, lrow_, lcol_, theta);
    IrecB = IFWTmatrix_m(IcompB, lrow_, lcol_, theta);

    IoutR_ = IrecR(1:Nor, 1:Mor);
    IoutG_ = IrecG(1:Nor, 1:Mor);
    IoutB_ = IrecB(1:Nor, 1:Mor);

    if initial_transformation
        IoutR_ = ITR2D(IoutR_, theta);
        IoutG_ = ITR2D(IoutG_, theta);
        IoutB_ = ITR2D(IoutB_, theta);
    end

    DR = abs(double(R) - double(IoutR_)).^2;
    DG = abs(double(G) - double(IoutG_)).^2;
    DB = abs(double(B) - double(IoutB_)).^2;
    mseVAL  = sum(DR(:) + DG(:) + DB(:)) / (numel(R) + numel(G) + numel(B));
    ssimVAL = (ssim(R, IoutR_) + ssim(G, IoutG_) + ssim(B, IoutB_)) / 3;

end