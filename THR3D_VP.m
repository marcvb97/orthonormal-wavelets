function [thr, zero_el, Icomp] = THR3D_VP(Idec, Iscal, Iwav, keep, lrow, lcol, lslice, flag)
%----------------------------------------
% AIM: Threshold the wavelet coefficients of a 3D tensor decomposition.
% Tensor extension of THR (which operates on 2D matrix decompositions).
%
% INPUT:
% Idec   = cell array of decomposed tensors from DEC3D_VP
% Iscal  = vector of scaling (approximation) coefficients
% Iwav   = vector of all wavelet (detail) coefficients
% keep   = fraction of coefficients to retain (0 < keep <= 1)
% lrow   = number of decompositions on rows   (dim 2)
% lcol   = number of decompositions on columns (dim 1)
% lslice = number of decompositions on slices  (dim 3)
% flag   = 0: threshold all coefficients (scaling + wavelet)
%          1: threshold wavelet coefficients only (keep all scaling)
%
% OUTPUT:
% thr     = threshold value used
% zero_el = number of coefficients set to zero
% Icomp   = cell array (same structure as Idec) with thresholded coefficients
%----------------------------------------

kmax = max([lrow, lcol, lslice]);

%% Determine which coefficients enter the threshold calculation
if flag == 0
    all_coef = [Iscal; Iwav];   % threshold everything
else
    all_coef = Iwav;            % threshold wavelet only, keep all scaling
end

%% Find threshold so that fraction 'keep' of all_coef are retained
num_keep = round(keep * length(all_coef));
sorted   = sort(abs(all_coef), 'descend');
if num_keep >= 1
    thr = sorted(num_keep);
else
    thr = Inf;
end

%% Apply threshold level by level
Icomp   = Idec;
zero_el = 0;

for k = 1:kmax

    J = Idec{k};
    [N, M, P] = size(J);
    Jt = J;

    % Scaling block size at this level
    if k <= lcol,   n = N/3;  else,  n = N;  end
    if k <= lrow,   m = M/3;  else,  m = M;  end
    if k <= lslice, p = P/3;  else,  p = P;  end

    if k == kmax
        % Coarsest level: scaling and wavelet subbands treated separately

        % Scaling block
        scal = Jt(1:n, 1:m, 1:p);
        if flag == 0
            mask_s       = abs(scal) < thr;
            zero_el      = zero_el + sum(mask_s(:));
            scal(mask_s) = 0;
        end
        Jt(1:n, 1:m, 1:p) = scal;

        % Wavelet subbands: everything outside the scaling block
        wav_region                 = true(N, M, P);
        wav_region(1:n, 1:m, 1:p) = false;
        to_zero                    = (abs(Jt) < thr) & wav_region;
        zero_el                    = zero_el + sum(to_zero(:));
        Jt(to_zero)                = 0;

    else
        % Finer levels: threshold the entire tensor block.
        % The scaling sub-block here is NOT the final approximation —
        % it gets overwritten during reconstruction by the coarser level,
        % so it is correct to threshold it here too.
        to_zero     = abs(Jt) < thr;
        zero_el     = zero_el + sum(to_zero(:));
        Jt(to_zero) = 0;
    end

    Icomp{k} = Jt;
end

end
