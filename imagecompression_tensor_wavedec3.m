function [RESULTS_GLOBAL, Ifin] = imagecompression_tensor_wavedec3(filepath)

% AIM: Tensor (3D) image compression using MATLAB built-in wavedec3/waverec3.
% The best result is found by varying wname (see vwname).
% Depending on the value of flag, all decomposition coefficients
% (wavelet and scaling) are compressed or only the wavelet coefficients.
% Results are displayed for all possible decomposition steps (see Kmax).

%% Settings
flag   = 1;     % 0 = compress all coefficients, 1 = compress details only

vkeep  = [0.5, 0.25, 0.10, 0.05, 0.02, 0.01, 0.005];
vkeep  = [0.005];
vwname = {'haar', 'db2', 'db4', 'sym4'};
vwname = {'db2'};

%% Load and prepare tensor
% Accepts either:
%   - a .mat file containing a 3D variable named 'I' (N x M x P double)
%   - a folder of grayscale images (one per slice)
[~, ~, ext] = fileparts(filepath);
if strcmpi(ext, '.mat')
    S  = load(filepath);
    fn = fieldnames(S);
    I  = double(S.(fn{1}));   % use first variable found
else
    % Load a folder of grayscale frames
    files = dir(fullfile(filepath, '*.png'));
    if isempty(files), files = dir(fullfile(filepath, '*.bmp')); end
    if isempty(files), files = dir(fullfile(filepath, '*.jpg')); end
    tmp = imread(fullfile(filepath, files(1).name));
    if size(tmp, 3) > 1, tmp = rgb2gray(tmp); end
    I = zeros(size(tmp,1), size(tmp,2), numel(files));
    for s = 1:numel(files)
        frame = imread(fullfile(filepath, files(s).name));
        if size(frame, 3) > 1, frame = rgb2gray(frame); end
        I(:,:,s) = double(frame);
    end
end

[Nor, Mor, Por] = size(I);
fprintf('Tensor size: %d x %d x %d\n', Nor, Mor, Por);

%% Decomposition parameters
Kmax = floor(log(min([Nor, Mor, Por])) / log(2) + 1.0e-12);
Kmax = 5;

%% Main loop
RESULTS   = zeros((Kmax-5+1) * length(vkeep) * length(vwname), 6);
TIMES     = zeros((Kmax-5+1) * length(vkeep) * length(vwname), 6);
cont      = 0;
cont_time = 0;

Ifin      = I;
err_best  = Inf;
wname_fin = vwname{1};
kmax_fin  = Kmax;

for kstep = 5:Kmax
    for keep = vkeep
        cont      = cont + 1;
        cont_time = cont_time + 1;

        %% First wavelet: initialise with first vwname value
        wname = vwname{1};

        tic;
        W = wavedec3(I, kstep, wname);
        t1 = toc;
        
        tic;
        [Wcomp, zero_el] = THR3D_wavedec(W, keep, flag, Nor*Mor*Por);
        t2 = toc;

        tic;
        Irec = waverec3(Wcomp);
        Iout = Irec(1:Nor, 1:Mor, 1:Por);
        t3 = toc;

        err       = sum(abs(I(:) - Iout(:)).^2) / numel(I);
        Ifin      = Iout;
        wname_fin = wname;
        kmax_fin  = kstep;
        zero_el_fin = zero_el;
        TIMES(cont_time, :) = [kstep, keep, 0, t1, t2, t3];

        %% Search remaining wavelet names for minimum MSE
        for wi = 2:length(vwname)
            wname = vwname{wi};
            cont_time = cont_time + 1;

            tic;
            W = wavedec3(I, kstep, wname);
            t1 = toc;

            tic;
            [Wcomp, zero_el_t] = THR3D_wavedec(W, keep, flag, Nor*Mor*Por);
            t2 = toc;

            tic;
            Irec = waverec3(Wcomp);
            Iout = Irec(1:Nor, 1:Mor, 1:Por);
            t3 = toc;

            mseVAL = sum(abs(I(:) - Iout(:)).^2) / numel(I);

            if mseVAL < err
                err       = mseVAL;
                Ifin      = Iout;
                wname_fin = wname;
                kmax_fin  = kstep;
                zero_el_fin = zero_el_t;
            end
            TIMES(cont_time, :) = [kstep, keep, wi, t1, t2, t3];
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

        num_total    = numel(I);
        pct_retained = (num_total - zero_el_fin) * 100 / num_total;

        RESULTS(cont, :) = [keep, psnrVAL, ssimVAL, pct_retained, kmax_fin, 0]

        % Display middle slice
        figure(cont)
        mid = round(Por / 2);
        imshow(uint8(Ifin(:,:,mid)))
        title(sprintf('Middle slice (%d/%d)  keep=%.3f  PSNR=%.2f dB  wavelet=%s', ...
            mid, Por, keep, psnrVAL, wname_fin))
        axis off equal
    end
end

%% Trim unused rows
RESULTS = RESULTS(1:cont, :);
TIMES   = TIMES(1:cont_time, :);

%% Display results table
varNames = {'CR', 'PSNR', 'SSIM', '% retained', 'Decomp levels', 'unused'};
disp('***** TENSOR COMPRESSION (wavedec3) *****')
T = array2table(RESULTS, 'VariableNames', varNames);
disp(T)

%% RESULTS_GLOBAL
RESULTS_GLOBAL         = struct;
RESULTS_GLOBAL.name    = filepath;
RESULTS_GLOBAL.RESULTS = RESULTS;
RESULTS_GLOBAL.TIMES   = TIMES;
RESULTS_GLOBAL.size    = [Nor, Mor, Por];

end


%% -------------------------------------------------------------------------
function [Wcomp, zero_el] = THR3D_wavedec(W, keep, flag, total_el)
% Threshold wavelet coefficients from wavedec3 output.
%
% W    : struct returned by wavedec3 (contains W.dec cell array)
% keep : fraction of coefficients to retain  (0 < keep <= 1)
% flag : 0 = threshold all coefficients
%        1 = threshold detail coefficients only (preserve approximation)
% total_el : total number of voxels, elements in the tensor
%
% Returns:
%   Wcomp    : thresholded struct (same format, ready for waverec3)
%   zero_el  : number of zeroed coefficients

Wcomp = W;

%% Identify approximation vs detail subbands
% wavedec3 stores coefficients in W.dec cell array
% Each cell contains one subband (7 detail + 1 approx per level)
% The approximation subband is W.dec{end}

n_subbands = length(W.dec);

if flag == 1
    % Compress detail coefficients only — leave approximation untouched
    detail_idx = 2:n_subbands;
else
    % Compress all coefficients including approximation
    detail_idx = 1:n_subbands;
end

% Gather all selected coefficients into one vector to find threshold
all_coefs = [];
for k = detail_idx
    all_coefs = [all_coefs; W.dec{k}(:)];
end

% Find threshold so that fraction 'keep' of coefficients are retained
% num_keep = round(keep * numel(all_coefs));
num_keep = round(keep * total_el - numel(W.dec{1}(:)));
sorted   = sort(abs(all_coefs), 'descend');
if num_keep >= 1
    thr = sorted(num_keep);
else
    thr = Inf;
end
thr

% Apply threshold
zero_el = 0;
for k = detail_idx
    d               = Wcomp.dec{k};
    mask            = abs(d) < thr;
    zero_el         = zero_el + sum(mask(:));
    d(mask)         = 0;
    Wcomp.dec{k}    = d;
end

end
