function Arec = IFWTmatrix3D_m(Adec, lrow, lcol, lslice, theta)
% IFWTmatrix3D_m  Multi-level 3D inverse wavelet transform.
%
% Reconstructs the tensor at full resolution from the cell array of
% scaling and wavelet coefficients produced by FWTmatrix3D_m.
%
% INPUT:
%   Adec   - cell array; Adec{k} contains transform coefficients at step k
%   lrow   - number of decompositions applied to rows   (dim 2)
%   lcol   - number of decompositions applied to columns (dim 1)
%   lslice - number of decompositions applied to slices  (dim 3)
%   theta  - localization parameter for the wavelet
%
% OUTPUT:
%   Arec   - reconstructed tensor (same size as Adec{1})

kmax = max([lrow, lcol, lslice]);
A    = Adec{kmax};

%% Inverse reconstruction loop (coarsest to finest)
for k = kmax:-1:1
    B         = Adec{k};
    [N, M, P] = size(B);

    if k <= lcol,   n = N / 3;  else,  n = N;  end
    if k <= lrow,   m = M / 3;  else,  m = M;  end
    if k <= lslice, p = P / 3;  else,  p = P;  end

    % Overwrite scaling block with reconstructed approximation
    B(1:n, 1:m, 1:p) = A(1:n, 1:m, 1:p);

    %-- Inverse slice transform (along dim 3) --%
    At = B;
    if k <= lslice
        for i = 1:N
            for j = 1:M
                % FT1step_m always outputs column vectors -> pass columns
                scaling = B(i, j, 1:p);    scaling = scaling(:);
                wavelet = B(i, j, p+1:P);  wavelet = wavelet(:);
                rec = IFT1step_m(scaling, wavelet, theta);  % 9x1
                At(i, j, 1:P) = rec(:);
            end
        end
    end

    %-- Inverse row transform (along dim 2) --%
    At2 = At;
    if k <= lrow
        for i = 1:N
            for p_idx = 1:P
                % Stored as row in Adec (transposed after FT1step_m)
                % -> read as row, convert to column for IFT1step_m
                scaling = At(i, 1:m,   p_idx);  scaling = scaling(:);
                wavelet = At(i, m+1:M, p_idx);  wavelet = wavelet(:);
                rec = IFT1step_m(scaling, wavelet, theta);  % 9x1
                At2(i, 1:M, p_idx) = rec.';     % store back as row
            end
        end
        At = At2;
    end

    %-- Inverse column transform (along dim 1) --%
    At3 = At;
    if k <= lcol
        for j = 1:M
            for p_idx = 1:P
                % Stored as column in Adec (FT1step_m output directly)
                % -> already a column, no transpose needed
                scaling = At(1:n,   j, p_idx);  scaling = scaling(:);
                wavelet = At(n+1:N, j, p_idx);  wavelet = wavelet(:);
                rec = IFT1step_m(scaling, wavelet, theta);  % 9x1
                At3(1:N, j, p_idx) = rec;        % store back as column
            end
        end
        At = At3;
    end

    A = At;
end

Arec = A;
end