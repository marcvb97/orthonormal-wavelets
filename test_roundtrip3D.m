%% test_roundtrip3D.m
% Round-trip test: DEC3D followed by IFWTmatrix3D_m should reconstruct
% the original tensor to within numerical precision.

clear; clc;
clear functions  % flush any cached old versions

%% Parameters
theta  = 0.5;   % localization parameter
kstep  = 4;     % maximum decomposition steps
dmin   = 3;     % minimum dimension before stopping

%% Test cases: various tensor sizes (including non-multiples of 3)
sizes = {
    [27,  27,  27],   % perfect powers of 3
    [24,  30,  18],   % multiples of 3
    [25,  31,  20],   % non-multiples of 3 (require padding)
    [10,  15,  12],   % small tensor
    [50,  40,  36],   % mixed sizes
};

fprintf('%-25s  %-12s  %-15s  %s\n', 'Tensor size', 'Max error', 'Rel error', 'Pass?');
fprintf('%s\n', repmat('-', 1, 65));

all_passed = true;

for t = 1:numel(sizes)
    sz = sizes{t};
    N  = sz(1);  M = sz(2);  P = sz(3);

    % Random input tensor
    rng(t);
    I = rand(N, M, P);

    %-- Forward transform --%
    [Idec, ~, ~, lrow, lcol, lslice] = DEC3D(I, theta, kstep, dmin);

    %-- Inverse transform --%
    Irec = IFWTmatrix3D_m(Idec, lrow, lcol, lslice, theta);

    % Trim back to original size (DEC3D may have padded)
    Irec = Irec(1:N, 1:M, 1:P);

    %-- Error metrics --%
    abs_err = max(abs(I(:) - Irec(:)));
    rel_err = abs_err / (max(abs(I(:))) + eps);
    passed  = abs_err < 1e-10;

    if ~passed
        all_passed = false;
    end

    fprintf('[%3d x %3d x %3d]%-4s  %-12.2e  %-15.2e  %s\n', ...
        N, M, P, '', abs_err, rel_err, yesno(passed));
end

fprintf('%s\n', repmat('-', 1, 65));
if all_passed
    fprintf('All tests PASSED.\n');
else
    fprintf('Some tests FAILED.\n');
end

%% Visual spot-check on a small example
fprintf('\n--- Visual spot-check (4x4x4 tensor) ---\n');
rng(42);
I_small = rand(4, 4, 4);
[Idec, ~, ~, lrow, lcol, lslice] = DEC3D(I_small, theta, kstep, dmin);
Irec_small = IFWTmatrix3D_m(Idec, lrow, lcol, lslice, theta);
Irec_small = Irec_small(1:4, 1:4, 1:4);

fprintf('Original slice 1:\n');
disp(I_small(:,:,1));
fprintf('Reconstructed slice 1:\n');
disp(Irec_small(:,:,1));
fprintf('Max absolute error: %.2e\n', max(abs(I_small(:) - Irec_small(:))));

%% Helper
function s = yesno(flag)
    if flag, s = 'PASS'; else, s = 'FAIL'; end
end