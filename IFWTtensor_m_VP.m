function Arec = IFWTtensor_m_VP(Adec, lrow, lcol, lslice, theta, fsca, fwav)
%----------------------------------------
% AIM: Reconstruct the approximation at the original resolution level (size)
% from the scaling and wavelet coefficients obtained by the 3D decompositions.
% This is the tensor (3D) extension of IFWTmatrix_m_VP.
%------------------------------------------
% INPUT:
% Adec   = cell array of transformed tensors
%          Adec{k} = full tensor of scal.+wav. coeff. at the k-th decomp. step
%                    size is always [Nj x Mj x Pj] with Nj,Mj,Pj multiples of 3
% lrow   = number of decompositions performed on rows   (dim 2)
% lcol   = number of decompositions performed on columns (dim 1)
% lslice = number of decompositions performed on slices  (dim 3)
% theta  = localization parameter for the wavelet
% fsca, fwav = normalization factors for the 1-step transforms
%
% OUTPUT:
% Arec = reconstructed tensor, same size as Adec{1}
%
% NB: The size of Adec{k} along each axis must always be a multiple of 3.
%     At each step k, if the reconstructed scaling tensor is larger than
%     needed it is trimmed (by 1 or 2 entries) to match the required size.
%----------------------------------------

kmax = max([lrow, lcol, lslice]);
A    = Adec{kmax};   % scaling+wavelet coefficients from the last decomposition

trow   = theta;
tcol   = theta;
tslice = theta;

%% Reconstruction loop (from coarsest to finest)
for k = kmax:-1:1

    %----------------------------------------------------------------------
    % INITIAL SETTING
    %----------------------------------------------------------------------
    B = Adec{k};
    [N, M, P] = size(B);

    % Size of the scaling part along each dimension at step k
    if k <= lcol
        n = N / 3;
    else
        n = N;   % columns were not decomposed at this step
    end

    if k <= lrow
        m = M / 3;
    else
        m = M;   % rows were not decomposed at this step
    end

    if k <= lslice
        p = P / 3;
    else
        p = P;   % slices were not decomposed at this step
    end

    % Insert current scaling coefficients into B
    B(1:n, 1:m, 1:p) = A(1:n, 1:m, 1:p);

    %----------------------------------------------------------------------
    % RECONSTRUCTION ORDER: reverse of forward transform
    % Forward was:  cols (dim1) → rows (dim2) → slices (dim3)
    % Inverse is:  slices (dim3) → rows (dim2) → cols (dim1)
    %----------------------------------------------------------------------

    %--- STEP 1: IFWT along SLICES (dim 3)
    At1 = B;
    if k <= lslice
        tslice = theta;
        for i = 1:N
            for j = 1:M
                scal = B(i, j, 1:p);    scal = scal(:);
                wav  = B(i, j, p+1:P);  wav  = wav(:);
                [rec, ~] = IFWT1step_m_VP(scal, wav, tslice, fsca, fwav);
                At1(i, j, 1:P) = rec;
            end
        end
    end

    %--- STEP 2: IFWT along ROWS (dim 2)
    At2 = At1;
    if k <= lrow
        trow = theta;
        for p_idx = 1:P
            for i = 1:N
                scal = At1(i, 1:m,   p_idx);  scal = scal(:);
                wav  = At1(i, m+1:M, p_idx);  wav  = wav(:);
                [rec, ~] = IFWT1step_m_VP(scal, wav, trow, fsca, fwav);
                At2(i, 1:M, p_idx) = rec.';
            end
        end
    end

    %--- STEP 3: IFWT along COLUMNS (dim 1)
    At3 = At2;
    if k <= lcol
        tcol = theta;
        for p_idx = 1:P
            for j = 1:M
                scal = At2(1:n,   j, p_idx);  scal = scal(:);
                wav  = At2(n+1:N, j, p_idx);  wav  = wav(:);
                [rec, ~] = IFWT1step_m_VP(scal, wav, tcol, fsca, fwav);
                At3(1:N, j, p_idx) = rec;
            end
        end
    end

    %----------------------------------------------------------------------
    % Final setting for next (finer) reconstruction step
    %----------------------------------------------------------------------
    A = At3;

end

%% Output
Arec = A;

end
