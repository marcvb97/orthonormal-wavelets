%% test_imagecompression_tensor_VP.m
% AIM: Test the full tensor compression pipeline:
%      DEC3D_VP -> THR3D_VP -> IFWTtensor_m_VP
%
% Tests cover:
%   1. Lossless case  (keep=1.0): reconstruction must equal original
%   2. Compression quality vs keep rate on a synthetic tensor
%   3. Compression quality vs theta on a synthetic tensor
%   4. Compression on a video-like tensor (uint8 range, 0-255)
%   5. Varying tensor sizes (padding paths)
%   6. flag=0 (threshold all) vs flag=1 (threshold details only)

%% Shared parameters
fsca  = sqrt(3);
fwav  = sqrt(3/2);
dmin  = 3;
tol   = 1e-10;   % for perfect reconstruction check

fprintf('=============================================================\n');
fprintf('  TENSOR COMPRESSION PIPELINE TEST\n');
fprintf('  DEC3D_VP -> THR3D_VP -> IFWTtensor_m_VP\n');
fprintf('=============================================================\n\n');

%% --------------------------------------------------------
%  Shared helper: run one compression and return metrics
%% --------------------------------------------------------
function [psnr_val, ssim_val, pct_ret, zero_el] = compress_tensor(I, theta, kstep, keep, dmin, fsca, fwav, flag)
    [Nor, Mor, Por] = size(I);

    [Idec, Iscal, Iwav, lrow, lcol, lslice] = DEC3D_VP(I, theta, kstep, dmin, fsca, fwav);
    [~, zero_el, Icomp]                      = THR3D_VP(Idec, Iscal, Iwav, keep, lrow, lcol, lslice, flag);
    Irec                                      = IFWTtensor_m_VP(Icomp, lrow, lcol, lslice, theta, fsca, fwav);
    Iout                                      = Irec(1:Nor, 1:Mor, 1:Por);

    err      = sum(abs(I(:) - Iout(:)).^2) / numel(I);
    if err > 0
        psnr_val = 10 * log10(255^2 / err);
    else
        psnr_val = Inf;
    end

    ssim_val = 0;
    for s = 1:Por
        ssim_val = ssim_val + ssim(I(:,:,s), Iout(:,:,s));
    end
    ssim_val = ssim_val / Por;

    num_total = length(Iscal) + length(Iwav);
    pct_ret   = (num_total - zero_el) * 100 / numel(I);
end

%% ============================================================
%  TEST 1: Perfect reconstruction (keep = 1.0)
%  All coefficients retained -> output must equal input exactly
%% ============================================================
fprintf('--- TEST 1: Perfect reconstruction (keep=1.0) ---\n');
theta = 0.5;  kstep = 3;  keep = 1.0;  flag = 1;
sizes = {[9,9,9], [27,27,9], [10,11,13]};

for si = 1:length(sizes)
    sz = sizes{si};
    A  = rand(sz(1), sz(2), sz(3)) * 255;
    [Nor, Mor, Por] = size(A);

    [Idec, Iscal, Iwav, lrow, lcol, lslice] = DEC3D_VP(A, theta, kstep, dmin, fsca, fwav);
    [~, ~, Icomp]                            = THR3D_VP(Idec, Iscal, Iwav, keep, lrow, lcol, lslice, flag);
    Irec                                      = IFWTtensor_m_VP(Icomp, lrow, lcol, lslice, theta, fsca, fwav);
    Iout                                      = Irec(1:Nor, 1:Mor, 1:Por);

    err = max(abs(A(:) - Iout(:)));
    if err < tol
        status = 'PASS';
    else
        status = 'FAIL';
    end
    fprintf('    Size %dx%dx%d  max_err=%.3e  --> %s\n', sz(1), sz(2), sz(3), err, status);
end
fprintf('\n');

%% ============================================================
%  TEST 2: Compression quality vs keep rate
%  Expect: lower keep -> lower PSNR, fewer retained coefficients
%% ============================================================
fprintf('--- TEST 2: Quality vs keep rate (27x27x9, theta=0.5, kstep=3, flag=1) ---\n');
fprintf('    %8s  %10s  %10s  %10s\n', 'keep', 'PSNR(dB)', 'SSIM', '% retained');

A     = rand(27, 27, 9) * 255;
theta = 0.5;  kstep = 3;  flag = 1;
vkeep = [1.0, 0.5, 0.25, 0.10, 0.05, 0.01];

prev_psnr = Inf;
all_ok    = true;
for keep = vkeep
    [psnr_val, ssim_val, pct_ret, ~] = compress_tensor(A, theta, kstep, keep, dmin, fsca, fwav, flag);
    fprintf('    %8.3f  %10.2f  %10.4f  %10.2f\n', keep, psnr_val, ssim_val, pct_ret);
    if psnr_val > prev_psnr + 1e-6   % PSNR should not increase as keep drops
        all_ok = false;
    end
    prev_psnr = psnr_val;
end
if all_ok
    fprintf('    --> PASS (PSNR decreases monotonically with keep)\n\n');
else
    fprintf('    --> FAIL (PSNR did not decrease monotonically)\n\n');
end

%% ============================================================
%  TEST 3: Compression quality vs theta
%  Shows which theta gives best PSNR for a fixed keep
%% ============================================================
fprintf('--- TEST 3: Quality vs theta (27x27x9, keep=0.1, kstep=3, flag=1) ---\n');
fprintf('    %8s  %10s  %10s\n', 'theta', 'PSNR(dB)', 'SSIM');

A      = rand(27, 27, 9) * 255;
keep   = 0.10;  kstep = 3;  flag = 1;
vtheta = 0.1:0.1:0.9;

best_psnr  = -Inf;
best_theta = NaN;
for theta = vtheta
    [psnr_val, ssim_val, ~, ~] = compress_tensor(A, theta, kstep, keep, dmin, fsca, fwav, flag);
    fprintf('    %8.2f  %10.2f  %10.4f\n', theta, psnr_val, ssim_val);
    if psnr_val > best_psnr
        best_psnr  = psnr_val;
        best_theta = theta;
    end
end
fprintf('    Best theta=%.2f  PSNR=%.2f dB\n\n', best_theta, best_psnr);

%% ============================================================
%  TEST 4: Video-like tensor (uint8 range 0-255)
%  Simulates a real grayscale video
%% ============================================================
fprintf('--- TEST 4: Video-like tensor (27x27x9, uint8 range) ---\n');
fprintf('    %8s  %10s  %10s  %10s\n', 'keep', 'PSNR(dB)', 'SSIM', '% retained');

A     = double(uint8(rand(27, 27, 9) * 255));
theta = 0.5;  kstep = 3;  flag = 1;
vkeep = [0.5, 0.25, 0.10, 0.05];

for keep = vkeep
    [psnr_val, ssim_val, pct_ret, ~] = compress_tensor(A, theta, kstep, keep, dmin, fsca, fwav, flag);
    fprintf('    %8.3f  %10.2f  %10.4f  %10.2f\n', keep, psnr_val, ssim_val, pct_ret);
end
fprintf('\n');

%% ============================================================
%  TEST 5: Various tensor sizes (padding paths)
%  Verifies the pipeline works for different dimension combinations
%% ============================================================
fprintf('--- TEST 5: Various tensor sizes ---\n');
fprintf('    %15s  %8s  %10s  %10s\n', 'Size', 'keep', 'PSNR(dB)', 'SSIM');

theta = 0.5;  kstep = 3;  keep = 0.10;  flag = 1;
sizes = {
    [9,  9,  9 ],   % all multiples of 3 - no padding
    [10, 10, 10],   % none multiple of 3 - all padded
    [27, 9,  3 ],   % mixed
    [9,  27, 6 ],   % mixed
    [12, 15, 10],   % all need padding
};

for si = 1:length(sizes)
    sz = sizes{si};
    A  = rand(sz(1), sz(2), sz(3)) * 255;
    [psnr_val, ssim_val, ~, ~] = compress_tensor(A, theta, kstep, keep, dmin, fsca, fwav, flag);
    fprintf('    %4dx%4dx%4d  %8.3f  %10.2f  %10.4f\n', sz(1), sz(2), sz(3), keep, psnr_val, ssim_val);
end
fprintf('\n');

%% ============================================================
%  TEST 6: flag=0 vs flag=1
%  flag=0: threshold all (including scaling)
%  flag=1: threshold details only (keep all scaling)
%  Expect: flag=1 gives better PSNR at same keep
%% ============================================================
fprintf('--- TEST 6: flag=0 vs flag=1 (27x27x9, theta=0.5, kstep=3) ---\n');
fprintf('    %8s  %10s  %10s  %10s\n', 'keep', 'PSNR flag0', 'PSNR flag1', 'flag1 better?');

A     = rand(27, 27, 9) * 255;
theta = 0.5;  kstep = 3;
vkeep = [0.25, 0.10, 0.05, 0.01];

all_ok = true;
for keep = vkeep
    [p0, ~, ~, ~] = compress_tensor(A, theta, kstep, keep, dmin, fsca, fwav, 0);
    [p1, ~, ~, ~] = compress_tensor(A, theta, kstep, keep, dmin, fsca, fwav, 1);
    better = p1 >= p0;
    if ~better, all_ok = false; end
    fprintf('    %8.3f  %10.2f  %10.2f  %10s\n', keep, p0, p1, mat2str(better));
end
if all_ok
    fprintf('    --> PASS (flag=1 always >= flag=0)\n\n');
else
    fprintf('    --> FAIL (flag=0 outperformed flag=1 in some case)\n\n');
end

%% ============================================================
%  TEST 7: Decomposition level comparison
%  Vary kstep for fixed keep and theta
%% ============================================================
fprintf('--- TEST 7: Quality vs decomposition level kstep (27x27x9, theta=0.5, keep=0.1) ---\n');
fprintf('    %8s  %10s  %10s\n', 'kstep', 'PSNR(dB)', 'SSIM');

A     = rand(27, 27, 9) * 255;
theta = 0.5;  keep = 0.10;  flag = 1;
Kmax  = floor(log(min(size(A))) / log(3) + 1e-12);

for kstep = 1:Kmax
    [psnr_val, ssim_val, ~, ~] = compress_tensor(A, theta, kstep, keep, dmin, fsca, fwav, flag);
    fprintf('    %8d  %10.2f  %10.4f\n', kstep, psnr_val, ssim_val);
end
fprintf('\n');

fprintf('=============================================================\n');
fprintf('  All tests complete.\n');
fprintf('=============================================================\n');
