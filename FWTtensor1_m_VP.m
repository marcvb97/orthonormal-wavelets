function [Adec, lrow, lcol, lslice] = FWTtensor1_m_VP(A, theta, kmax, dmin, fsca, fwav)
%----------------------------------------
% AIM: Decompose the 3D tensor A kmax times (or less) while the dimensions
% (n. of rows, n. of cols, n. of slices) of the scaling part to further
% decompose is > dmin.
% This is the tensor (3D) extension of FWTmatrix1_m_VP.
%------------------------------------------
% INPUT:
% A     = 3D tensor to decompose  [N1 x M1 x P1]
% theta = localization degree of the wavelet
% kmax  = max number of decompositions on rows, cols AND slices
% dmin  = min dimension along each axis to apply one more decomposition
% fsca, fwav = normalization factors for the 1-step transforms
%
% OUTPUT:
% Adec        = cell array: Adec{k} = full tensor of scal.+wav. coeff.
%                           at the k-th decomposition step  [Nj x Mj x Pj]
% lrow        = number of decompositions performed on rows
% lcol        = number of decompositions performed on columns
% lslice      = number of decompositions performed on slices
%
% STRUCTURE OF ONE DECOMPOSITION STEP (analogous to 2D):
%   Along columns (dim 1): FWT1step on each column-fiber  → nj scaling rows
%   Along rows    (dim 2): FWT1step on each row-fiber     → mj scaling cols
%   Along slices  (dim 3): FWT1step on each slice-fiber   → pj scaling slices
% scaling coefficients live in A(1:nj, 1:mj, 1:pj)
%----------------------------------------

%% Initialization
[N1, M1, P1] = size(A);
dmax  = max([N1, M1, P1]);
lrow  = 0;  lcol = 0;  lslice = 0;
kCONT = 0;

%% Decomposition loop
while dmax > dmin && kCONT < kmax
    kCONT = kCONT + 1;

    %----------------------------------------------------------------------
    % STEP 1: pad ROWS (dim 1) to Nj = 3*nj if needed
    %----------------------------------------------------------------------
    if N1 > dmin
        r1 = rem(N1, 3);
        if r1 == 0
            Nj = N1;  nj = N1/3;
        else
            nplus = 3 - r1;
            Nj    = N1 + nplus;
            nj    = Nj / 3;
            for k = 1:nplus
                A(N1+k, :, :) = A(N1, :, :);   % replicate last row
            end
        end
    else
        Nj = N1;  nj = N1;   % skip column-transform
    end

    %----------------------------------------------------------------------
    % STEP 2: pad COLUMNS (dim 2) to Mj = 3*mj if needed
    %----------------------------------------------------------------------
    if M1 > dmin
        r1 = rem(M1, 3);
        if r1 == 0
            Mj = M1;  mj = M1/3;
        else
            nplus = 3 - r1;
            Mj    = M1 + nplus;
            mj    = Mj / 3;
            for k = 1:nplus
                A(:, M1+k, :) = A(:, M1, :);   % replicate last column
            end
        end
    else
        Mj = M1;  mj = M1;   % skip row-transform
    end

    %----------------------------------------------------------------------
    % STEP 3: pad SLICES (dim 3) to Pj = 3*pj if needed
    %----------------------------------------------------------------------
    if P1 > dmin
        r1 = rem(P1, 3);
        if r1 == 0
            Pj = P1;  pj = P1/3;
        else
            nplus = 3 - r1;
            Pj    = P1 + nplus;
            pj    = Pj / 3;
            for k = 1:nplus
                A(:, :, P1+k) = A(:, :, P1);   % replicate last slice
            end
        end
    else
        Pj = P1;  pj = P1;   % skip slice-transform
    end

    %----------------------------------------------------------------------
    % STEP 4: FWT along COLUMNS (dim 1) — transform each fiber A(:, j, p)
    %----------------------------------------------------------------------
    At1 = A;   % will hold result after column transform
    if Nj > dmin
        lcol = lcol + 1;
        tcol = theta;
        for p = 1:Pj
            for j = 1:Mj
                fiber = A(1:Nj, j, p);  fiber = fiber(:);   % force to column vector
                [anew, bnew, ~] = FWT1step_m_VP(fiber, tcol, fsca, fwav);
                At1(1:Nj, j, p) = [anew; bnew];
            end
        end
    end

    %----------------------------------------------------------------------
    % STEP 5: FWT along ROWS (dim 2) — transform each fiber At1(i, :, p)
    %----------------------------------------------------------------------
    At2 = At1;   % will hold result after row transform
    if Mj > dmin
        lrow = lrow + 1;
        trow = theta;
        for p = 1:Pj
            for i = 1:Nj
                fiber = At1(i, 1:Mj, p);  fiber = fiber(:);   % force to column vector
                [anew, bnew, ~] = FWT1step_m_VP(fiber, trow, fsca, fwav);
                At2(i, 1:Mj, p) = [anew; bnew].';
            end
        end
    end

    %----------------------------------------------------------------------
    % STEP 6: FWT along SLICES (dim 3) — transform each fiber At2(i, j, :)
    %----------------------------------------------------------------------
    At3 = At2;   % will hold result after slice transform
    if Pj > dmin
        lslice = lslice + 1;
        tslice = theta;
        for i = 1:Nj
            for j = 1:Mj
                fiber = At2(i, j, 1:Pj);  fiber = fiber(:);   % force to column vector
                [anew, bnew, ~] = FWT1step_m_VP(fiber, tslice, fsca, fwav);
                At3(i, j, 1:Pj) = [anew; bnew];
            end
        end
    end

    %----------------------------------------------------------------------
    % STEP 7: store result and prepare for next decomposition
    %----------------------------------------------------------------------
    Adec{kCONT} = At3(1:Nj, 1:Mj, 1:Pj);   % full coeff tensor at step k

    % Scaling part only → to be decomposed next iteration
    A  = At3(1:nj, 1:mj, 1:pj);
    N1 = nj;  M1 = mj;  P1 = pj;
    dmax = max([N1, M1, P1]);

end
end
