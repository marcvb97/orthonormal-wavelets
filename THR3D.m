function [thres, zero_el, Icomp] = THR3D(Idec, Iscal, Iwav, keep, lrow, lcol, lslice, flag)
% THR3D  Threshold-based compression of 3D wavelet coefficients.
%
% INPUT:
%   Idec    - cell array of transform tensors (from FWTmatrix3D_m)
%   Iscal   - vector of scaling (approximation) coefficients
%   Iwav    - vector of all wavelet (detail) coefficients
%   keep    - fraction of total coefficients to retain (0 < keep <= 1)
%   lrow    - number of row decompositions
%   lcol    - number of column decompositions
%   lslice  - number of slice decompositions
%   flag    - 0: threshold all coefficients (scaling + wavelet)
%             1: threshold wavelet coefficients only (preserve scaling)
%
% OUTPUT:
%   thres   - threshold value applied
%   zero_el - number of coefficients set to zero
%   Icomp   - cell array of compressed transform tensors

kmax = max([lrow, lcol, lslice]);

if flag == 0
    %% Compress all coefficients (scaling + wavelet)
    Itot    = [Iscal; Iwav];
    stot    = sort(abs(Itot), 'descend');
    num     = numel(Itot);
    ii      = floor(num * keep);
    thres   = stot(ii);
    zero_el = numel(find(Itot .* (abs(Itot) < thres)));

    for k = kmax:-1:1
        A        = Idec{k};
        Icomp{k} = A .* (abs(A) >= thres);
    end

else
    %% Compress wavelet coefficients only (scaling part is always kept)
    swav = sort(abs(Iwav), 'descend');
    Nsca = length(Iscal);
    Nwav = length(Iwav);
    num  = Nsca + Nwav;
    ii   = floor(num * keep - Nsca);

    if ii > Nwav || ii < 1
        % Nothing to threshold: retain everything
        thres   = 1e3;
        zero_el = 0;
        Icomp   = Idec;
    else
        thres   = swav(ii);
        zero_el = numel(find(Iwav .* (abs(Iwav) < thres)));

        for k = kmax:-1:1
            A         = Idec{k};
            [N, M, P] = size(A);

            % Threshold all coefficients
            It = A .* (abs(A) > thres);

            % Determine size of scaling block at this level
            % matching the same logic as DEC3D/IFWTmatrix3D_m
            if k <= lcol
                n = N / 3;
            else
                n = N;
            end

            if k <= lrow
                m = M / 3;
            else
                m = M;
            end

            if k <= lslice
                p = P / 3;
            else
                p = P;
            end

            % Restore the scaling block (never thresholded)
            It(1:n, 1:m, 1:p) = A(1:n, 1:m, 1:p);
            Icomp{k}           = It;
        end
    end
end
end