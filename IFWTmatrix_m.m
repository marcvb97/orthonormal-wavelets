function Arec = IFWTmatrix_m(Adec, lrow, lcol, theta)
% IFWTmatrix_m  Multi-level 2D inverse wavelet transform.
%
% Reconstructs the image at full resolution from the cell array of
% scaling and wavelet coefficients produced by FWTmatrix1_m.
%
% INPUT:
%   Adec   - cell array; Adec{k} contains transform coefficients at step k
%   lrow   - number of decompositions applied to rows
%   lcol   - number of decompositions applied to columns
%   theta  - localization parameter for the wavelet
%
% OUTPUT:
%   Arec   - reconstructed matrix (same size as Adec{1})
%
% NOTE:
%   Each Adec{k} must have size [3*n, 3*m].  If the reconstructed
%   scaling matrix has more than n rows (or m columns) after a step,
%   the excess rows (columns) are discarded before the next step.

kmax = max(lrow, lcol);
A    = Adec{kmax};

%% Inverse reconstruction loop (coarsest to finest)
for k = kmax:-1:1

    B      = Adec{k};
    [N, M] = size(B);

    % Size of the scaling submatrix at this level
    if k <= lcol
        n = N/3;
    else
        n = N;      % columns were not decomposed at this level
    end
    if k <= lrow
        m = M/3;
    else
        m = M;      % rows were not decomposed at this level
    end

    % Overwrite the scaling block with the reconstructed approximation
    B(1:n, 1:m) = A(1:n, 1:m);

    %-- Inverse row transform --%
    if k <= lrow
        trow = theta;
        for i = 1:N
            [At(i, 1:M), ~] = IFT1step_m(B(i, 1:m).', B(i, m+1:M).', trow);
        end
    else
        At = B;
    end

    %-- Inverse column transform --%
    if k <= lcol
        tcol = theta;
        for j = 1:M
            [At(1:N, j), ~] = IFT1step_m(At(1:n, j), At(n+1:N, j), tcol);
        end
    end

    A = At;
end

Arec = A;

end
