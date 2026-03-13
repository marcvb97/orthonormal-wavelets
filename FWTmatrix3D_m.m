function [Adec, lrow, lcol, lslice] = FWTmatrix3D_m(A, theta, kmax, dmin)
% FWTmatrix3D_m  Multi-level 3D forward wavelet transform.
%
% Decomposes 3D tensor A up to kmax times, stopping early when the scaling
% sub-tensor dimension falls to dmin or below.
%
% INPUT:
%   A      - 3D tensor to decompose (N1 x M1 x P1)
%   theta  - localization parameter for the wavelet
%   kmax   - maximum number of decomposition steps
%   dmin   - minimum dimension of the scaling part before stopping
%
% OUTPUT:
%   Adec   - cell array; Adec{k} contains all scaling and wavelet
%            coefficients produced at step k (size always [3*n, 3*m, 3*p])
%   lrow   - number of decompositions applied to rows   (dim 2)
%   lcol   - number of decompositions applied to columns (dim 1)
%   lslice - number of decompositions applied to slices  (dim 3)

[N1, M1, P1] = size(A);
dmax   = max([N1, M1, P1]);
lrow   = 0;
lcol   = 0;
lslice = 0;
kCONT  = 0;

while dmax > dmin && kCONT < kmax
    kCONT = kCONT + 1;

    %-- Pad dim 1 (rows) to a multiple of 3 --%
    if N1 > dmin
        r1 = rem(N1, 3);
        if r1 == 0
            Nj = N1;
        else
            nplus = 3 - r1;
            Nj = N1 + nplus;
            for k = 1:nplus
                A(N1+k, :, :) = A(N1, :, :);
            end
        end
        nj = Nj / 3;
    else
        Nj = N1;  nj = N1;
    end

    %-- Pad dim 2 (columns) to a multiple of 3 --%
    if M1 > dmin
        r1 = rem(M1, 3);
        if r1 == 0
            Mj = M1;
        else
            nplus = 3 - r1;
            Mj = M1 + nplus;
            for k = 1:nplus
                A(:, M1+k, :) = A(:, M1, :);
            end
        end
        mj = Mj / 3;
    else
        Mj = M1;  mj = M1;
    end

    %-- Pad dim 3 (slices) to a multiple of 3 --%
    if P1 > dmin
        r1 = rem(P1, 3);
        if r1 == 0
            Pj = P1;
        else
            nplus = 3 - r1;
            Pj = P1 + nplus;
            for k = 1:nplus
                A(:, :, P1+k) = A(:, :, P1);
            end
        end
        pj = Pj / 3;
    else
        Pj = P1;  pj = P1;
    end

    %-- Transform along dim 1 (columns): read A, write At1 --%
    At1 = A;
    if Nj > dmin
        lcol = lcol + 1;
        for p = 1:Pj
            for j = 1:Mj
                [anew, bnew, ~] = FT1step_m(A(1:Nj, j, p), theta);  % read A
                At1(1:Nj, j, p) = [anew; bnew];
            end
        end
    end

    %-- Transform along dim 2 (rows): read At1, write At2 --%
    At2 = At1;
    if Mj > dmin
        lrow = lrow + 1;
        for p = 1:Pj
            for i = 1:Nj
                row = At1(i, 1:Mj, p);
                [anew, bnew, ~] = FT1step_m(row(:), theta);  % (:) forces column
                At2(i, 1:Mj, p) = [anew; bnew].';
            end
        end
    end

    %-- Transform along dim 3 (slices): read At2, write At3 --%
    At3 = At2;
    if Pj > dmin
        lslice = lslice + 1;
        for i = 1:Nj
            for j = 1:Mj
                slice = At2(i, j, 1:Pj);
                [anew, bnew, ~] = FT1step_m(slice(:), theta);  % (:) forces column
                At3(i, j, 1:Pj) = [anew; bnew].';
            end
        end
    end
    %-- Prepare for next iteration --%
    A           = At3(1:nj, 1:mj, 1:pj);
    Adec{kCONT} = At3(1:Nj, 1:Mj, 1:Pj);
    N1 = nj;  M1 = mj;  P1 = pj;
    dmax = max([N1, M1, P1]);
end
end