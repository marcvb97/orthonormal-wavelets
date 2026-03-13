function [Idec, Iscal, Iwav, lrow, lcol, lslice] = DEC3D(I, theta, kstep, dmin)
% DEC3D  Wavelet decomposition of 3D tensor I.
%
% Pads I to have dimensions that are multiples of 3, applies the 3D forward
% wavelet transform, then extracts scaling and wavelet coefficient vectors.
%
% INPUT:
%   I      - input 3D tensor (N0 x M0 x P0)
%   theta  - localization parameter for the wavelet
%   kstep  - maximum number of decomposition steps
%   dmin   - minimum dimension before stopping decomposition
%
% OUTPUT:
%   Idec   - cell array of transform tensors at each decomposition step
%   Iscal  - vector of scaling (approximation) coefficients
%   Iwav   - vector of all wavelet (detail) coefficients
%   lrow   - number of decompositions applied to rows
%   lcol   - number of decompositions applied to columns
%   lslice - number of decompositions applied to slices

[N0, M0, P0] = size(I);
Iin = I;

%% Pad rows to a multiple of 3
r0 = rem(N0, 3);
if r0 ~= 0
    Rnplus = 3 - r0;
    for k = 1:Rnplus
        Iin(N0+k, :, :) = Iin(N0, :, :);
    end
end

%% Pad columns to a multiple of 3
r0 = rem(M0, 3);
if r0 ~= 0
    Cnplus = 3 - r0;
    for k = 1:Cnplus
        Iin(:, M0+k, :) = Iin(:, M0, :);
    end
end

%% Pad slices to a multiple of 3
r0 = rem(P0, 3);
if r0 ~= 0
    Snplus = 3 - r0;
    for k = 1:Snplus
        Iin(:, :, P0+k) = Iin(:, :, P0);
    end
end

%% Apply 3D forward wavelet transform
[Idec, lrow, lcol, lslice] = FWTmatrix3D_m(Iin, theta, kstep, dmin);

%% Extract scaling and wavelet coefficient vectors from the coarsest level
kmax = max([lrow, lcol, lslice]);
J    = Idec{kmax};
[N, M, P] = size(J);

% Determine row/col/slice scaling subblock sizes
if kmax <= lrow
    n = N / 3;
else
    n = N;
end

if kmax <= lcol
    m = M / 3;
else
    m = M;
end

if kmax <= lslice
    p = P / 3;
else
    p = P;
end

% Scaling coefficients: the (n x m x p) corner block
Jscal = J(1:n, 1:m, 1:p);
Iscal = Jscal(:);

% Wavelet coefficients: everything outside the scaling block
Iwav = extract_wav3d(J, n, m, p, N, M, P);

%% Collect wavelet coefficients from all coarser decomposition levels
if kmax > 1
    for k = kmax-1:-1:1
        J         = Idec{k};
        [N, M, P] = size(J);

        if k <= lrow
            n = N / 3;
        else
            n = N;
        end

        if k <= lcol
            m = M / 3;
        else
            m = M;
        end

        if k <= lslice
            p = P / 3;
        else
            p = P;
        end

        Iwav = [Iwav; extract_wav3d(J, n, m, p, N, M, P)];
    end
end

end % function


%% Helper: extract all wavelet coefficients from a 3D coefficient tensor
function Iwav = extract_wav3d(J, n, m, p, N, M, P)
% Collects all 7 wavelet subbands, i.e. everything outside the
% scaling block J(1:n, 1:m, 1:p). The 8 subbands of a 3D wavelet
% step are (S=scaling, W=wavelet along that axis):
%
%   SSS : J(1:n,   1:m,   1:p  )  <- scaling block, excluded
%   SSW : J(1:n,   1:m,   p+1:P)
%   SWS : J(1:n,   m+1:M, 1:p  )
%   SWW : J(1:n,   m+1:M, p+1:P)
%   WSS : J(n+1:N, 1:m,   1:p  )
%   WSW : J(n+1:N, 1:m,   p+1:P)
%   WWS : J(n+1:N, m+1:M, 1:p  )
%   WWW : J(n+1:N, m+1:M, p+1:P)

SSW = J(1:n,   1:m,   p+1:P);
SWS = J(1:n,   m+1:M, 1:p  );
SWW = J(1:n,   m+1:M, p+1:P);
WSS = J(n+1:N, 1:m,   1:p  );
WSW = J(n+1:N, 1:m,   p+1:P);
WWS = J(n+1:N, m+1:M, 1:p  );
WWW = J(n+1:N, m+1:M, p+1:P);

Iwav = [SSW(:); SWS(:); SWW(:); WSS(:); WSW(:); WWS(:); WWW(:)];

% Handle degenerate cases where some axes were not transformed
% (the corresponding wavelet subband will be empty and is safely ignored)
end