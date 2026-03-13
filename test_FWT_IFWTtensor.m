%% test_FWT_IFWTtensor.m
% AIM: Test the perfect reconstruction property of FWTtensor1_m_VP and
%      IFWTtensor_m_VP on several synthetic tensors.
%
% A perfect wavelet transform satisfies:
%      IFWTtensor( FWTtensor( A ) ) == A   (up to floating point precision)
%
% We test:
%   1. Small tensor (all dimensions multiple of 3)        → no padding needed
%   2. Small tensor (dimensions NOT multiples of 3)       → padding triggered
%   3. Asymmetric tensor (different sizes per dimension)  → mixed padding
%   4. Larger random tensor (stress test)
%   5. Constant tensor                                    → trivial signal
%   6. Single-slice tensor (P=3, quasi-2D)
%
% For each test we report:
%   - max absolute error  (should be < 1e-10 for perfect reconstruction)
%   - lrow, lcol, lslice  (number of decompositions per axis)
%   - PASS / FAIL

%% Parameters shared across all tests
theta = 0.5;
fsca  = sqrt(3);
fwav  = sqrt(3/2);
kmax  = 4;
dmin  = 3;
tol   = 1e-10;   % tolerance for perfect reconstruction

fprintf('=============================================================\n');
fprintf('  PERFECT RECONSTRUCTION TEST: FWTtensor + IFWTtensor\n');
fprintf('  theta=%.2f  fsca=%.4f  fwav=%.4f  kmax=%d  dmin=%d\n', ...
        theta, fsca, fwav, kmax, dmin);
fprintf('=============================================================\n\n');

all_pass = true;

%% --------------------------------------------------------
%  Helper: run one test and print result
%% --------------------------------------------------------
function run_test(label, A, theta, kmax, dmin, fsca, fwav, tol)
    fprintf('--- %s ---\n', label);
    fprintf('    Input size : %d x %d x %d\n', size(A,1), size(A,2), size(A,3));

    % Forward transform
    [Adec, lrow, lcol, lslice] = FWTtensor1_m_VP(A, theta, kmax, dmin, fsca, fwav);
    fprintf('    lrow=%d  lcol=%d  lslice=%d\n', lrow, lcol, lslice);
    fprintf('    Decomp levels: %d\n', length(Adec));

    % Inverse transform
    Arec = IFWTtensor_m_VP(Adec, lrow, lcol, lslice, theta, fsca, fwav);

    % Trim reconstructed tensor to original size (padding may have enlarged it)
    [N, M, P] = size(A);
    Arec_trim = Arec(1:N, 1:M, 1:P);

    % Error
    err = max(abs(A(:) - Arec_trim(:)));
    fprintf('    Max abs error: %.3e\n', err);

    if err < tol
        fprintf('    --> PASS\n\n');
    else
        fprintf('    --> FAIL  (error exceeds tolerance %.1e)\n\n', tol);
    end
end

%% --------------------------------------------------------
%  TEST 1: Small tensor, all dims multiple of 3
%% --------------------------------------------------------
A1 = rand(9, 9, 9);
run_test('TEST 1: 9x9x9 random (all dims mult of 3)', ...
         A1, theta, kmax, dmin, fsca, fwav, tol);

%% --------------------------------------------------------
%  TEST 2: Small tensor, dims NOT multiples of 3
%% --------------------------------------------------------
A2 = rand(10, 11, 13);
run_test('TEST 2: 10x11x13 random (dims NOT mult of 3)', ...
         A2, theta, kmax, dmin, fsca, fwav, tol);

%% --------------------------------------------------------
%  TEST 3: Asymmetric tensor
%% --------------------------------------------------------
A3 = rand(9, 27, 6);
run_test('TEST 3: 9x27x6 asymmetric random', ...
         A3, theta, kmax, dmin, fsca, fwav, tol);

%% --------------------------------------------------------
%  TEST 4: Larger stress test
%% --------------------------------------------------------
A4 = rand(30, 40, 20);
run_test('TEST 4: 30x40x20 larger random tensor', ...
         A4, theta, kmax, dmin, fsca, fwav, tol);

%% --------------------------------------------------------
%  TEST 5: Constant tensor (trivial signal)
%% --------------------------------------------------------
A5 = 128 * ones(9, 9, 9);
run_test('TEST 5: 9x9x9 constant tensor (all 128)', ...
         A5, theta, kmax, dmin, fsca, fwav, tol);

%% --------------------------------------------------------
%  TEST 6: Quasi-2D tensor (single group of slices, P=3)
%% --------------------------------------------------------
A6 = rand(9, 9, 3);
run_test('TEST 6: 9x9x3 quasi-2D (P=3)', ...
         A6, theta, kmax, dmin, fsca, fwav, tol);

%% --------------------------------------------------------
%  TEST 7: Multiple theta values
%% --------------------------------------------------------
fprintf('--- TEST 7: varying theta on 9x9x9 ---\n');
A7 = rand(9, 9, 9);
for theta_t = [0.1, 0.3, 0.5, 0.7, 0.9]
    [Adec7, lrow7, lcol7, lslice7] = FWTtensor1_m_VP(A7, theta_t, kmax, dmin, fsca, fwav);
    Arec7      = IFWTtensor_m_VP(Adec7, lrow7, lcol7, lslice7, theta_t, fsca, fwav);
    Arec7_trim = Arec7(1:9, 1:9, 1:9);
    err7       = max(abs(A7(:) - Arec7_trim(:)));
    if err7 < tol
        status = 'PASS';
    else
        status = 'FAIL';
    end
    fprintf('    theta=%.1f  max_err=%.3e  --> %s\n', theta_t, err7, status);
end
fprintf('\n');

%% --------------------------------------------------------
%  TEST 8: Energy conservation check
%  The wavelet transform should preserve signal energy
%% --------------------------------------------------------
fprintf('--- TEST 8: energy conservation on 9x9x9 ---\n');
A8 = rand(9, 9, 9);
[Adec8, lrow8, lcol8, lslice8] = FWTtensor1_m_VP(A8, theta, kmax, dmin, fsca, fwav);

% Collect all wavelet coefficients from the last decomposition level
coef_last = Adec8{end}(:);
energy_in  = sum(A8(:).^2);
energy_dec = sum(coef_last.^2);
ratio      = energy_dec / energy_in;
fprintf('    Energy input      : %.6f\n', energy_in);
fprintf('    Energy (last dec) : %.6f\n', energy_dec);
fprintf('    Ratio             : %.6f  (close to 1 if near-isometric)\n\n', ratio);

%% --------------------------------------------------------
%  TEST 9: Video-like tensor (grayscale frames 0-255)
%% --------------------------------------------------------
A9 = double(uint8(rand(27, 27, 9) * 255));
run_test('TEST 9: 27x27x9 uint8-range video-like tensor', ...
         A9, theta, kmax, dmin, fsca, fwav, tol);

fprintf('=============================================================\n');
fprintf('  All tests complete.\n');
fprintf('=============================================================\n');
