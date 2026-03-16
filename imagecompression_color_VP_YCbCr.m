function [RESULTS_GLOBAL,Ifin2] = imagecompression_color_VP(filepath)

% AIM: Image compression using custom wavelets.
% The best result is found by varying theta (see vtheta).
% Depending on the value of flag, all decomposition coefficients 
% (wavelet and scaling) are compressed or only the wavelet coefficients.
% Results are displayed for all possible decomposition steps (see Kmax).
% initial_transformation indicates if TR2D and ITR2D are used.

%% Settings
flag = 1;   % 0 = compress all coefficients, 1 = compress details only
threshold_type = 1; % 1 = separately, 2 = combined over the channels

vkeep  = [0.5, 0.25, 0.10, 0.05, 0.02, 0.01, 0.005];
% vkeep  = [0.005];
vtheta = [0.1:0.1:0.9];
vtheta = [0.5];

% normalization factors for the 1step transforms
% NB in decomposing we multiply the output by these factors
% in reconstructing we divide the input by these factors
fsca = sqrt(3); fwav = sqrt(3/2); 

%% Load and prepare image
% testImage = '../PEXELS300/997704.bmp';
% testImage = '../13US/fig01.png';
% testImage = '../KODAK/1.png';
testImage = filepath;
II = imread(testImage);
imshow(II);
ycbcr_img = rgb2ycbcr(II);
Y  = ycbcr_img(:,:,1); % Luma (intensity)
Cb = ycbcr_img(:,:,2); % Blue chroma
Cr = ycbcr_img(:,:,3); % Red chroma
% we are denoting Y by R, Cb by G, and Cr by B
R = double(Y); G = double(Cb); B = double(Cr);

[Nor, Mor] = size(R);

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
        cont = cont + 1;
	cont_time = cont_time+1;

        %% First theta: initialise error with first vtheta value
        theta = vtheta(1);

        AR = R;
        AG = G;
        AB = B;
	tic;
        [IdecR, IscalR, IwavR, lrow, lcol] = DEC_VP(AR, theta, kstep, dmin, fsca, fwav);
        [IdecG, IscalG, IwavG, lrow, lcol] = DEC_VP(AG, theta, kstep, dmin, fsca, fwav);
        [IdecB, IscalB, IwavB, lrow, lcol] = DEC_VP(AB, theta, kstep, dmin, fsca, fwav);
	t1 = toc;
	tic;
        switch threshold_type
            case 1
                coeff_kept = keep*3*Nor*Mor;
                coeff_kept_Y = 0.60*coeff_kept;
                coeff_kept_C = 0.20*coeff_kept;
                keep_Y = coeff_kept_Y / (Nor*Mor);
                keep_C = coeff_kept_C / (Nor*Mor);
                [thr, zero_elR, IcompR] = THR(IdecR, IscalR, IwavR, keep_Y, lrow, lcol, flag);
                [thr, zero_elG, IcompG] = THR(IdecG, IscalG, IwavG, keep_C, lrow, lcol, flag);
                [thr, zero_elB, IcompB] = THR(IdecB, IscalB, IwavB, keep_C, lrow, lcol, flag);
                zero_el = zero_elR + zero_elG + zero_elB;
            case 2
                [thr, zero_el, IcompR, IcompG, IcompB] = THR_color(IdecR, IdecG, IdecB,...
                    IscalR, IscalG, IscalB, IwavR, IwavG, IwavB, keep, ...
                    lrow, lcol, flag);
        end
	t2 = toc;
	tic;
        IrecR = IFWTmatrix_m_VP(IcompR, lrow, lcol, theta, fsca, fwav);
        IrecG = IFWTmatrix_m_VP(IcompG, lrow, lcol, theta, fsca, fwav);
        IrecB = IFWTmatrix_m_VP(IcompB, lrow, lcol, theta, fsca, fwav);

        IoutR = IrecR(1:Nor, 1:Mor);
        IoutG = IrecG(1:Nor, 1:Mor);
        IoutB = IrecB(1:Nor, 1:Mor);
	t3 = toc;
        DR    = abs(double(R) - double(IoutR)).^2;
        DG    = abs(double(G) - double(IoutG)).^2;
        DB    = abs(double(B) - double(IoutB)).^2;
        err  = sum(DR(:)+DG(:)+DB(:)) / (numel(R)+numel(G)+numel(B));

        Rfin      = IoutR;
        Gfin      = IoutG;
        Bfin      = IoutB;
        Ifin      = cat(3, Rfin, Gfin, Bfin);
        theta_fin = theta;
        kmax      = max(lrow, lcol);
        TIMES(cont_time,:) = [kstep, keep, theta, t1, t2, t3, err];


        %% Search remaining theta values for minimum MSE
        for theta = setdiff(vtheta, vtheta(1))
            cont_time = cont_time+1;
            AR = R;
            AG = G;
            AB = B;
	    tic;
            [IdecR, IscalR, IwavR, lrow, lcol] = DEC_VP(AR, theta, kstep, dmin, fsca, fwav);
            [IdecG, IscalG, IwavG, lrow, lcol] = DEC_VP(AG, theta, kstep, dmin, fsca, fwav);
            [IdecB, IscalB, IwavB, lrow, lcol] = DEC_VP(AB, theta, kstep, dmin, fsca, fwav);
            t1 = toc;
            tic;
            switch threshold_type
                case 1
                    coeff_kept = keep*3*Nor*Mor;
                    coeff_kept_Y = 0.60*coeff_kept;
                    coeff_kept_C = 0.20*coeff_kept;
                    keep_Y = coeff_kept_Y / (Nor*Mor);
                    keep_C = coeff_kept_C / (Nor*Mor);
                    [thr, zero_elR, IcompR] = THR(IdecR, IscalR, IwavR, keep_Y, lrow, lcol, flag);
                    [thr, zero_elG, IcompG] = THR(IdecG, IscalG, IwavG, keep_C, lrow, lcol, flag);
                    [thr, zero_elB, IcompB] = THR(IdecB, IscalB, IwavB, keep_C, lrow, lcol, flag);
                    zero_el = zero_elR + zero_elG + zero_elB;
                case 2
                    [thr, zero_el, IcompR, IcompG, IcompB] = THR_color(IdecR, IdecG, IdecB,...
                        IscalR, IscalG, IscalB, IwavR, IwavG, IwavB, keep, ...
                        lrow, lcol, flag);
            end
	    t2 = toc;
            tic;
            IrecR = IFWTmatrix_m_VP(IcompR, lrow, lcol, theta, fsca, fwav);
            IrecG = IFWTmatrix_m_VP(IcompG, lrow, lcol, theta, fsca, fwav);
            IrecB = IFWTmatrix_m_VP(IcompB, lrow, lcol, theta, fsca, fwav);

            IoutR = IrecR(1:Nor, 1:Mor);
            IoutG = IrecG(1:Nor, 1:Mor);
            IoutB = IrecB(1:Nor, 1:Mor);
            t3 = toc;
            DR    = abs(double(R) - double(IoutR)).^2;
            DG    = abs(double(G) - double(IoutG)).^2;
            DB    = abs(double(B) - double(IoutB)).^2;
            mseVAL  = sum(DR(:)+DG(:)+DB(:)) / (numel(R)+numel(G)+numel(B));

            if mseVAL < err
                err       = mseVAL;
                Rfin      = IoutR;
                Gfin      = IoutG;
                Bfin      = IoutB;
                Ifin      = cat(3, Rfin, Gfin, Bfin);
                theta_fin = theta;
                kmax      = max(lrow, lcol);
            end
            TIMES(cont_time,:) = [kstep, keep, theta, t1, t2, t3, mseVAL];
        end

        %% Quality metrics
        disp(filepath)
        psnrVAL = 10 * log10(255^2 / err);
        SSIMVAL  = (ssim(R, Rfin)+ssim(G, Gfin)+ssim(B, Bfin))/3;

        num = 3*length([IscalR; IwavR]);

        RESULTS(cont, :) = [keep, psnrVAL, SSIMVAL, (num - zero_el)*100/(3*Nor*Mor), theta_fin, kmax]

        figure(cont)
        Ifin2 = uint8(ycbcr2rgb(uint8(Ifin)));
        imshow(Ifin2);
        axis off equal
    end
end

%% Display results table
RESULTS_best = zeros(length(vkeep),size(RESULTS,2));
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
RESULTS_GLOBAL.name = filepath;
RESULTS_GLOBAL.RESULTS = RESULTS;
RESULTS_GLOBAL.RESULTS_best = RESULTS_best;
RESULTS_GLOBAL.TIMES = TIMES;
RESULTS_GLOBAL.size = [Nor, Mor];

end