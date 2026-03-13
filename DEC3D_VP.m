function [Idec, Iscal, Iwav, lrow, lcol, lslice] = DEC3D_VP(I, theta, kstep, dmin, fsca, fwav)
%----------------------------------------
% AIM: Decompose the 3D tensor I and separate the result into:
%   - Iscal: vector of scaling (approximation) coefficients
%   - Iwav:  vector of all wavelet (detail) coefficients
% This is the tensor (3D) extension of DEC_VP.
%
% INPUT:
% I      = 3D tensor [N0 x M0 x P0]
% theta  = localization parameter for the wavelet
% kstep  = max number of decomposition steps
% dmin   = min dimension to allow a further decomposition
% fsca, fwav = normalization factors
%
% OUTPUT:
% Idec   = cell array of decomposed tensors (from FWTtensor1_m_VP)
% Iscal  = vector of approximation coefficients (from last decomp. level)
% Iwav   = vector of ALL wavelet coefficients (all levels, coarsest first)
% lrow   = number of decompositions on rows   (dim 2)
% lcol   = number of decompositions on columns (dim 1)
% lslice = number of decompositions on slices  (dim 3)
%
% SUBBAND STRUCTURE per decomposition step (analogous to 2D):
% In 2D: 1 scaling (LL) + 3 wavelet subbands (LH, HL, HH)
% In 3D: 1 scaling (LLL) + 7 wavelet subbands:
%        LLH, LHL, LHH, HLL, HLH, HHL, HHH
% where L=low(scaling) H=high(wavelet) along dims [col, row, slice]
% Scaling part lives in J(1:n, 1:m, 1:p)
%----------------------------------------

[N0, M0, P0] = size(I);
Iin = I;

%% Pad each dimension to a multiple of 3
% pad rows (dim 1)
r0 = rem(N0, 3);
if r0 ~= 0
    Rnplus = 3 - r0;
    for k = 1:Rnplus
        Iin(N0+k, :, :) = Iin(N0, :, :);
    end
end

% pad columns (dim 2)
r0 = rem(M0, 3);
if r0 ~= 0
    Cnplus = 3 - r0;
    for k = 1:Cnplus
        Iin(:, M0+k, :) = Iin(:, M0, :);
    end
end

% pad slices (dim 3)
r0 = rem(P0, 3);
if r0 ~= 0
    Snplus = 3 - r0;
    for k = 1:Snplus
        Iin(:, :, P0+k) = Iin(:, :, P0);
    end
end

%% Apply 3D decomposition
[Idec, lrow, lcol, lslice] = FWTtensor1_m_VP(Iin, theta, kstep, dmin, fsca, fwav);

%% Separate scaling and wavelet coefficients
% At each level k, the full tensor Idec{k} has size [N x M x P] where
% N, M, P are multiples of 3. The scaling part is J(1:n, 1:m, 1:p).
% The 7 wavelet subbands are the remaining 7 block regions.

kmax = max([lrow, lcol, lslice]);   % actual number of decompositions done

%% Helper: extract the 7 wavelet subbands from tensor J at level k
% n,m,p = size of scaling part along each dim
% Returns them concatenated as a column vector

    function wav_k = extract_wav_subbands(J, k, lcol, lrow, lslice)
        [N, M, P] = size(J);

        % Determine scaling size along each dim
        if k <= lcol,   n = N/3;  else,  n = N;  end
        if k <= lrow,   m = M/3;  else,  m = M;  end
        if k <= lslice, p = P/3;  else,  p = P;  end

        % Index ranges
        rL_N = 1:n;      rH_N = n+1:N;   % col  (dim1): L=scaling H=wavelet
        rL_M = 1:m;      rH_M = m+1:M;   % row  (dim2)
        rL_P = 1:p;      rH_P = p+1:P;   % slice(dim3)

        wav_k = [];

        % Only include subbands for dimensions that were actually decomposed
        col_dec   = (k <= lcol);
        row_dec   = (k <= lrow);
        slice_dec = (k <= lslice);

        if col_dec && row_dec && slice_dec
            % All 7 non-LLL subbands
            LLH = J(rL_N, rL_M, rH_P);   % col=L row=L slice=H
            LHL = J(rL_N, rH_M, rL_P);   % col=L row=H slice=L
            LHH = J(rL_N, rH_M, rH_P);   % col=L row=H slice=H
            HLL = J(rH_N, rL_M, rL_P);   % col=H row=L slice=L
            HLH = J(rH_N, rL_M, rH_P);   % col=H row=L slice=H
            HHL = J(rH_N, rH_M, rL_P);   % col=H row=H slice=L
            HHH = J(rH_N, rH_M, rH_P);   % col=H row=H slice=H
            wav_k = [LLH(:); LHL(:); LHH(:); HLL(:); HLH(:); HHL(:); HHH(:)];

        elseif col_dec && row_dec && ~slice_dec
            % 2D-like: 3 subbands (LH, HL, HH) — no slice wavelet
            LH = J(rL_N, rH_M, rL_P);
            HL = J(rH_N, rL_M, rL_P);
            HH = J(rH_N, rH_M, rL_P);
            wav_k = [LH(:); HL(:); HH(:)];

        elseif col_dec && ~row_dec && slice_dec
            % col + slice only
            LH_slice = J(rL_N, rL_M, rH_P);
            HL_slice = J(rH_N, rL_M, rL_P);
            HH_slice = J(rH_N, rL_M, rH_P);
            wav_k = [LH_slice(:); HL_slice(:); HH_slice(:)];

        elseif ~col_dec && row_dec && slice_dec
            % row + slice only
            LH_slice = J(rL_N, rL_M, rH_P);
            HL_slice = J(rL_N, rH_M, rL_P);
            HH_slice = J(rL_N, rH_M, rH_P);
            wav_k = [LH_slice(:); HL_slice(:); HH_slice(:)];

        elseif col_dec && ~row_dec && ~slice_dec
            % col only
            wav_k = J(rH_N, rL_M, rL_P);
            wav_k = wav_k(:);

        elseif ~col_dec && row_dec && ~slice_dec
            % row only
            wav_k = J(rL_N, rH_M, rL_P);
            wav_k = wav_k(:);

        elseif ~col_dec && ~row_dec && slice_dec
            % slice only
            wav_k = J(rL_N, rL_M, rH_P);
            wav_k = wav_k(:);
        end
    end

%% Coarsest level: extract Iscal and first Iwav
J = Idec{kmax};
[N, M, P] = size(J);

if kmax <= lcol,   n = N/3;  else,  n = N;  end
if kmax <= lrow,   m = M/3;  else,  m = M;  end
if kmax <= lslice, p = P/3;  else,  p = P;  end

Jscal = J(1:n, 1:m, 1:p);
Iscal = Jscal(:);
Iwav  = extract_wav_subbands(J, kmax, lcol, lrow, lslice);

%% Finer levels: accumulate wavelet coefficients
if kmax > 1
    for k = kmax-1:-1:1
        J     = Idec{k};
        Iwav  = [Iwav; extract_wav_subbands(J, k, lcol, lrow, lslice)];
    end
end

end
