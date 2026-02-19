function [Idec, Iscal, Iwav, lrow, lcol] = DEC(I, theta, kstep, dmin)
% DEC  Wavelet decomposition of image I.
%
% Pads I to have dimensions that are multiples of 3, applies the forward
% wavelet transform, then extracts scaling and wavelet coefficient vectors.
%
% INPUT:
%   I      - input image matrix
%   theta  - localization parameter for the wavelet
%   kstep  - maximum number of decomposition steps
%   dmin   - minimum dimension before stopping decomposition
%
% OUTPUT:
%   Idec   - cell array of transform matrices at each decomposition step
%   Iscal  - vector of scaling (approximation) coefficients
%   Iwav   - vector of all wavelet (detail) coefficients
%   lrow   - number of decompositions applied to rows
%   lcol   - number of decompositions applied to columns

[N0, M0] = size(I);
Iin = I;

%% Pad rows to a multiple of 3
r0 = rem(N0, 3);
if r0 ~= 0
    Rnplus = 3 - r0;
    for k = 1:Rnplus
        Iin(N0+k, :) = Iin(N0, :);
    end
end

%% Pad columns to a multiple of 3
r0 = rem(M0, 3);
if r0 ~= 0
    Cnplus = 3 - r0;
    for k = 1:Cnplus
        Iin(:, M0+k) = Iin(:, M0);
    end
end

%% Apply forward wavelet transform
[Idec, lrow, lcol] = FWTmatrix1_m(Iin, theta, kstep, dmin);

%% Extract scaling and wavelet coefficient vectors
kmax = max(lrow, lcol);
J    = Idec{kmax};
[N, M] = size(J);

if kmax <= lcol && kmax <= lrow
    n = N/3;  m = M/3;
    I1   = J(1:n,   m+1:M);
    I2   = J(n+1:N, 1:m);
    I3   = J(n+1:N, m+1:M);
    Iwav = [I1(:); I2(:); I3(:)];
elseif kmax <= lcol && kmax > lrow
    n    = N/3;  m = M;
    I1   = J(n+1:N, 1:M);
    Iwav = I1(:);
elseif kmax > lcol && kmax <= lrow
    n    = N;  m = M/3;
    I1   = J(1:n, m+1:M);
    Iwav = I1(:);
end

Jscal = J(1:n, 1:m);
Iscal = Jscal(:);

%% Collect wavelet coefficients from coarser decomposition levels
if kmax > 1
    for k = kmax-1:-1:1
        J      = Idec{k};
        [N, M] = size(J);

        if k <= lcol && k <= lrow
            n = N/3;  m = M/3;
            I1   = J(1:n,   m+1:M);
            I2   = J(n+1:N, 1:m);
            I3   = J(n+1:N, m+1:M);
            Iwav = [Iwav; I1(:); I2(:); I3(:)];
        elseif k <= lcol && k > lrow
            n    = N/3;  m = M;
            I1   = J(n+1:N, 1:M);
            Iwav = [Iwav; I1(:)];
        elseif k > lcol && k <= lrow
            n    = N;  m = M/3;
            I1   = J(1:n, m+1:M);
            Iwav = [Iwav; I1(:)];
        end
    end
end

end
