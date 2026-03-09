function [RESULTS_GLOBAL] = imagecompression_color_VP_par(filepath)

% AIM: Image compression using custom wavelets.
% The best result is found by varying theta (see vtheta).
% Depending on the value of flag, all decomposition coefficients 
% (wavelet and scaling) are compressed or only the wavelet coefficients.
% Results are displayed for all possible decomposition steps (see Kmax).
% initial_transformation indicates if TR2D and ITR2D are used.
% R, G, B channel computations are parallelized using parfor.

%% Settings
flag = 1;   % 0 = compress all coefficients, 1 = compress details only
threshold_type = 1; % 1 = separately, 2 = combined over the channels

vkeep  = [0.5, 0.25, 0.10, 0.05, 0.02, 0.01, 0.005];
% vkeep  = [0.5];
vtheta = [0.1:0.1:0.9];
vtheta = [0.5];

% normalization factors for the 1step transforms
% NB in decomposing we multiply the output by these factors
% in reconstructing we divide the input by these factors
fsca = sqrt(3); fwav = sqrt(3/2);

%% Load and prepare image
testImage = filepath;
II = imread(testImage);
imshow(II)
channels{1} = double(II(:,:,1));  % R
channels{2} = double(II(:,:,2));  % G
channels{3} = double(II(:,:,3));  % B

[Nor, Mor] = size(channels{1});

%% Decomposition parameters
Kmax = floor(log(min(Nor, Mor))/log(3) + 1.0e-12);
% Kmax = min(5,Kmax);
Kmax = 4;
dmin = 3;

%% Main loop
RESULTS = zeros((Kmax-4+1) * length(vkeep), 6);
TIMES   = zeros((Kmax-4+1) * length(vkeep) * length(vtheta), 7);
cont = 0;
cont_time = 0;

for kstep = 4:Kmax
    for keep = vkeep
        cont      = cont + 1;
        cont_time = cont_time + 1;

        %% First theta: initialise error with first vtheta value
        theta = vtheta(1);

        % --- Parallel decomposition + threshold + reconstruction per channel ---
        t1_arr   = zeros(1,3);
        t2_arr   = zeros(1,3);
        t3_arr   = zeros(1,3);
        Iout     = cell(1,3);
        zero_els = zeros(1,3);
        lrow_arr = zeros(1,3);
        lcol_arr = zeros(1,3);
        IscalLen = zeros(1,3);

        parfor ch = 1:3
            Chan = channels{ch};

            tic;
            [Idec, Iscal, Iwav, lrow, lcol] = DEC_VP(Chan, theta, kstep, dmin, fsca, fwav);
            t1_arr(ch) = toc;

            tic;
            [~, zero_el_ch, Icomp] = THR(Idec, Iscal, Iwav, keep, lrow, lcol, flag);
            t2_arr(ch) = toc;

            tic;
            Irec   = IFWTmatrix_m_VP(Icomp, lrow, lcol, theta, fsca, fwav);
            IoutCh = Irec(1:Nor, 1:Mor);
            t3_arr(ch) = toc;

            Iout{ch}     = IoutCh;
            zero_els(ch) = zero_el_ch;
            lrow_arr(ch) = lrow;
            lcol_arr(ch) = lcol;
            IscalLen(ch) = length([Iscal; Iwav]);
        end

        % threshold_type == 2 (combined) couples all channels; keep sequential.
        if threshold_type == 2
            [IdecR, IscalR, IwavR, lrow, lcol] = DEC_VP(channels{1}, theta, kstep, dmin, fsca, fwav);
            [IdecG, IscalG, IwavG, ~,    ~   ] = DEC_VP(channels{2}, theta, kstep, dmin, fsca, fwav);
            [IdecB, IscalB, IwavB, ~,    ~   ] = DEC_VP(channels{3}, theta, kstep, dmin, fsca, fwav);
            [~, zero_el, IcompR, IcompG, IcompB] = THR_color(IdecR, IdecG, IdecB, ...
                IscalR, IscalG, IscalB, IwavR, IwavG, IwavB, keep, lrow, lcol, flag);
            Iout{1} = IFWTmatrix_m_VP(IcompR, lrow, lcol, theta, fsca, fwav); Iout{1} = Iout{1}(1:Nor,1:Mor);
            Iout{2} = IFWTmatrix_m_VP(IcompG, lrow, lcol, theta, fsca, fwav); Iout{2} = Iout{2}(1:Nor,1:Mor);
            Iout{3} = IFWTmatrix_m_VP(IcompB, lrow, lcol, theta, fsca, fwav); Iout{3} = Iout{3}(1:Nor,1:Mor);
            zero_els = [zero_el/3, zero_el/3, zero_el/3];
        end

        t1 = max(t1_arr);  t2 = max(t2_arr);  t3 = max(t3_arr);
        zero_el = sum(zero_els);

        DR  = abs(double(channels{1}) - double(Iout{1})).^2;
        DG  = abs(double(channels{2}) - double(Iout{2})).^2;
        DB  = abs(double(channels{3}) - double(Iout{3})).^2;
        err = sum(DR(:)+DG(:)+DB(:)) / (numel(channels{1})*3);

        Ifin      = cat(3, Iout{1}, Iout{2}, Iout{3});
        theta_fin = theta;
        kmax      = max(lrow_arr(1), lcol_arr(1));
        TIMES(cont_time,:) = [kstep, keep, theta, t1, t2, t3, err];

        %% Search remaining theta values for minimum MSE
        for theta = setdiff(vtheta, vtheta(1))
            cont_time = cont_time + 1;

            t1_arr      = zeros(1,3);
            t2_arr      = zeros(1,3);
            t3_arr      = zeros(1,3);
            Iout_th     = cell(1,3);
            zero_els_th = zeros(1,3);

            parfor ch = 1:3
                Chan = channels{ch};

                tic;
                [Idec, Iscal, Iwav, lrow, lcol] = DEC_VP(Chan, theta, kstep, dmin, fsca, fwav);
                t1_arr(ch) = toc;

                tic;
                [~, zero_el_ch, Icomp] = THR(Idec, Iscal, Iwav, keep, lrow, lcol, flag);
                t2_arr(ch) = toc;

                tic;
                Irec   = IFWTmatrix_m_VP(Icomp, lrow, lcol, theta, fsca, fwav);
                IoutCh = Irec(1:Nor, 1:Mor);
                t3_arr(ch) = toc;

                Iout_th{ch}     = IoutCh;
                zero_els_th(ch) = zero_el_ch;
            end

            t1 = max(t1_arr);  t2 = max(t2_arr);  t3 = max(t3_arr);

            DR     = abs(double(channels{1}) - double(Iout_th{1})).^2;
            DG     = abs(double(channels{2}) - double(Iout_th{2})).^2;
            DB     = abs(double(channels{3}) - double(Iout_th{3})).^2;
            mseVAL = sum(DR(:)+DG(:)+DB(:)) / (numel(channels{1})*3);

            if mseVAL < err
                err       = mseVAL;
                Ifin      = cat(3, Iout_th{1}, Iout_th{2}, Iout_th{3});
                theta_fin = theta;
                kmax      = max(lrow_arr(1), lcol_arr(1));
                zero_el   = sum(zero_els_th);
                Iout      = Iout_th;
            end
            TIMES(cont_time,:) = [kstep, keep, theta, t1, t2, t3, mseVAL];
        end

        %% Quality metrics
        disp(filepath)
        psnrVAL = 10 * log10(255^2 / err);
        SSIMVAL  = (ssim(channels{1}, Iout{1}) + ...
                    ssim(channels{2}, Iout{2}) + ...
                    ssim(channels{3}, Iout{3})) / 3;

        num = 3 * IscalLen(1);  % same length for all channels

        RESULTS(cont, :) = [keep, psnrVAL, SSIMVAL, (num - zero_el)*100/(3*Nor*Mor), theta_fin, kmax]

        figure(cont)
        imshow(uint8(Ifin))
        axis off equal
    end
end

%% Display results table
RESULTS_best = zeros(length(vkeep), size(RESULTS,2));
for i = 1:length(vkeep)
    I = find(RESULTS(:,1) == vkeep(i));
    [~,J] = max(RESULTS(I,2));
    RESULTS_best(i,:) = RESULTS(I(J(1)),:);
end

varNames = {'CR', 'PSNR', 'SSIM', '\% retained', 'Best $\theta$', 'Decomp'};

disp('***** OUR CASE *****')
T = array2table(RESULTS_best, 'VariableNames', varNames);
disp(T)

%% RESULTS_GLOBAL
RESULTS_GLOBAL = struct;
RESULTS_GLOBAL.name         = filepath;
RESULTS_GLOBAL.RESULTS      = RESULTS;
RESULTS_GLOBAL.RESULTS_best = RESULTS_best;
RESULTS_GLOBAL.TIMES        = TIMES;
RESULTS_GLOBAL.size         = [Nor, Mor];

end
