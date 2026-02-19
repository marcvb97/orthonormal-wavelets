function [Adec, lrow, lcol] = FWTmatrix1_m(A, theta, kmax, dmin)
% FWTmatrix1_m  Multi-level 2D forward wavelet transform.
%
% Decomposes matrix A up to kmax times, stopping early when the scaling
% submatrix dimension falls to dmin or below.
%
% INPUT:
%   A      - matrix to decompose
%   theta  - localization parameter for the wavelet
%   kmax   - maximum number of decomposition steps
%   dmin   - minimum dimension of the scaling part before stopping
%
% OUTPUT:
%   Adec   - cell array; Adec{k} contains all scaling and wavelet
%            coefficients produced at step k (size always [3*n, 3*m])
%   lrow   - number of decompositions applied to rows
%   lcol   - number of decompositions applied to columns
%
% NOTE: Scaling coefficients occupy Adec{k}(1:n, 1:m);
%       wavelet coefficients occupy the remaining entries.

[N1, M1] = size(A);
dmax = max(N1, M1);

lrow  = 0;
lcol  = 0;
kCONT = 0;

%% Decomposition loop
while dmax > dmin && kCONT < kmax

    kCONT = kCONT + 1;

    %-- Pad rows to a multiple of 3 (if above dmin) --%
    if N1 > dmin
        n1 = floor(N1/3);
        r1 = rem(N1, 3);
        if r1 == 0
            Nj = N1;  nj = n1;
        else
            nplus = 3 - r1;
            Nj    = N1 + nplus;  nj = Nj/3;
            for k = 1:nplus
                A(N1+k, :) = A(N1, :);
            end
        end
    else
        Nj = N1;  nj = N1;   % columns too short: skip column transform
    end

    %-- Pad columns to a multiple of 3 (if above dmin) --%
    if M1 > dmin
        m1 = floor(M1/3);
        r1 = rem(M1, 3);
        if r1 == 0
            Mj = M1;  mj = m1;
        else
            nplus = 3 - r1;
            Mj    = M1 + nplus;  mj = Mj/3;
            for k = 1:nplus
                A(:, M1+k) = A(:, M1);
            end
        end
    else
        Mj = M1;  mj = M1;   % rows too short: skip row transform
    end

    %-- 1-step column transform --%
    if Nj > dmin
        lcol = lcol + 1;
        tcol = theta;
        for j = 1:Mj
            [anew, bnew, ~] = FT1step_m(A(1:Nj, j), tcol);
            At1(1:Nj, j)    = [anew; bnew];
        end
    else
        At1 = A;
    end

    %-- 1-step row transform --%
    if Mj > dmin
        lrow = lrow + 1;
        trow = theta;
        for i = 1:Nj
            [anew, bnew, ~] = FT1step_m(At1(i, 1:Mj).', trow);
            At(i, 1:Mj)     = [anew; bnew].';
        end
    else
        At = At1;
    end

    %-- Prepare for next iteration --%
    A          = At(1:nj, 1:mj);       % scaling submatrix to further decompose
    Adec{kCONT} = At(1:Nj, 1:Mj);    % full coefficient matrix at this step
    N1   = nj;  M1 = mj;
    dmax = max(N1, M1);
end

end
