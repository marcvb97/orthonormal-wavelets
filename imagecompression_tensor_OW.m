function [RESULTS_GLOBAL,Ifin] = imagecompression_tensor_OW(filepath)

% AIM: Tensor (3D) image compression using custom wavelets.
% The best result is found by varying theta (see vtheta).
% Depending on the value of flag, all decomposition coefficients
% (wavelet and scaling) are compressed or only the wavelet coefficients.
% Results are displayed for all possible decomposition steps (see Kmax).

%% Settings
flag               = 1;     % 0 = compress all coefficients, 1 = compress details only
initial_transformation = true;

vkeep  = [0.5, 0.25, 0.10, 0.05, 0.02, 0.01, 0.005];
vkeep  = [0.001];
vtheta = [0.1:0.1:0.9];
vtheta = [0.5];

%% Load and prepare tensor
% Accepts either:
%   - a .mat file containing a 3D variable named 'I' (N x M x P double)
%   - a folder of grayscale images (one per slice)
[~, ~, ext] = fileparts(filepath);
if strcmpi(ext, '.mat')
    S = load(filepath);
    fn = fieldnames(S);
    I = double(S.(fn{1}));   % use first variable found
else
    % Load a folder of grayscale frames
    files = dir(fullfile(filepath, '*.png'));
    if isempty(files), files = dir(fullfile(filepath, '*.bmp')); end
    if isempty(files), files = dir(fullfile(filepath, '*.jpg')); end
    tmp = imread(fullfile(filepath, files(1).name));
    if size(tmp,3) > 1, tmp = rgb2gray(tmp); end
    I = zeros(size(tmp,1), size(tmp,2), numel(files));
    for s = 1:numel(files)
        frame = imread(fullfile(filepath, files(s).name));
        if size(frame,3) > 1, frame = rgb2gray(frame); end
        I(:,:,s) = double(frame);
    end
end

[Nor, Mor, Por] = size(I);
fprintf('Tensor size: %d x %d x %d\n', Nor, Mor, Por);

%% Decomposition parameters
Kmax = floor(log(min([Nor, Mor, Por])) / log(3) + 1.0e-12);
Kmax = 4;
dmin = 3;

%% Main loop
RESULTS    = zeros((Kmax-3+1) * length(vkeep), 6);
TIMES      = zeros((Kmax-3+1) * length(vkeep) * length(vtheta), 7);
cont       = 0;
cont_time  = 0;

for kstep = 4:Kmax
    for keep = vkeep
        cont      = cont + 1;
        cont_time = cont_time + 1;

        %% First theta: initialise with first vtheta value
        theta = vtheta(1);

        tic;
        if initial_transformation
            AT = TR3D(I, theta);
        else
            AT = I;
        end
        [Idec, Iscal, Iwav, lrow, lcol, lslice] = DEC3D(AT, theta, kstep, dmin);
        t1 = toc;

        tic;
        [thr, zero_el, Icomp] = THR3D(Idec, Iscal, Iwav, keep, lrow, lcol, lslice, flag);
        t2 = toc;

        tic;
        Irec = IFWTmatrix3D_m(Icomp, lrow, lcol, lslice, theta);
        Iout = Irec(1:Nor, 1:Mor, 1:Por);
        if initial_transformation
            Iout = ITR3D(Iout, theta);
        end
        t3 = toc;

        err       = sum(abs(I(:) - Iout(:)).^2) / numel(I);
        Ifin      = Iout;
        theta_fin = theta;
        kmax_fin  = max([lrow, lcol, lslice]);
        TIMES(cont_time, :) = [kstep, keep, theta, t1, t2, t3, err];

        %% Search remaining theta values for minimum MSE
        for theta = setdiff(vtheta, vtheta(1))
            cont_time = cont_time + 1;

            tic;
            if initial_transformation
                AT = TR3D(I, theta);
            else
                AT = I;
            end
            [Idec, Iscal, Iwav, lrow, lcol, lslice] = DEC3D(AT, theta, kstep, dmin);
            t1 = toc;

            tic;
            [thr, zero_el_t, Icomp] = THR3D(Idec, Iscal, Iwav, keep, lrow, lcol, lslice, flag);
            t2 = toc;

            tic;
            Irec = IFWTmatrix3D_m(Icomp, lrow, lcol, lslice, theta);
            Iout = Irec(1:Nor, 1:Mor, 1:Por);
            if initial_transformation
                Iout = ITR3D(Iout, theta);
            end
            t3 = toc;

            mseVAL = sum(abs(I(:) - Iout(:)).^2) / numel(I);

            if mseVAL < err
                err       = mseVAL;
                Ifin      = Iout;
                theta_fin = theta;
                kmax_fin  = max([lrow, lcol, lslice]);
                zero_el   = zero_el_t;
            end
            TIMES(cont_time, :) = [kstep, keep, theta, t1, t2, t3, mseVAL];
        end

        %% Quality metrics
        disp(filepath)
        psnrVAL = 10 * log10(255^2 / err);

        % SSIM averaged over slices
        ssimVAL = 0;
        for s = 1:Por
            ssimVAL = ssimVAL + ssim(I(:,:,s), Ifin(:,:,s));
        end
        ssimVAL = ssimVAL / Por;

        num_total = numel(I);
        pct_retained = (num_total - zero_el) * 100 / num_total;

        RESULTS(cont, :) = [keep, psnrVAL, ssimVAL, pct_retained, theta_fin, kmax_fin]

        % Display middle slice
        figure(cont)
        mid = round(Por / 2);
        imshow(uint8(Ifin(:,:,mid)))
        title(sprintf('Middle slice (%d/%d)  keep=%.3f  PSNR=%.2f dB', mid, Por, keep, psnrVAL))
        axis off equal
    end
end

%% Display results table
RESULTS_best = zeros(length(vkeep), size(RESULTS, 2));
for i = 1:length(vkeep)
    I_idx = find(RESULTS(:,1) == vkeep(i));
    [~, J] = max(RESULTS(I_idx, 2));
    RESULTS_best(i, :) = RESULTS(I_idx(J(1)), :);
end

varNames = {'CR', 'PSNR', 'SSIM', '% retained', 'Best theta', 'Decomp'};
disp('***** TENSOR COMPRESSION *****')
T = array2table(RESULTS_best, 'VariableNames', varNames);
disp(T)

%% RESULTS_GLOBAL
RESULTS_GLOBAL = struct;
RESULTS_GLOBAL.name         = filepath;
RESULTS_GLOBAL.RESULTS      = RESULTS;
RESULTS_GLOBAL.RESULTS_best = RESULTS_best;
RESULTS_GLOBAL.TIMES        = TIMES;
RESULTS_GLOBAL.size         = [Nor, Mor, Por];

end