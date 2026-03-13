function [anew, bnew, nplus] = FWT1step_m_VP_fast(a, theta, fsca, fwav)
% INPUT:
%   a      - column vector to transform, length N1
%   theta  - parameter in ]0,n[ defining degree parameter m
%   fsca   - normalization factor for scaling coefficients
%   fwav   - normalization factor for wavelet coefficients
% OUTPUT:
%   anew   - scaling coefficients (column vector)
%   bnew   - wavelet coefficients (column vector)
%   nplus  - number of elements added to input vector a (0 if none)

    %------------------
    % Initialization
    %------------------
    N1 = length(a);
    r1 = rem(N1, 3);

    if r1 ~= 0
        nplus = 3 - r1;
        a(end+1:end+nplus) = a(N1);  % vectorized padding
    else
        nplus = 0;
    end

    N  = length(a);
    n  = N / 3;

    % Degree parameter m
    if theta < 1
        m = max(1, floor(theta * n));   % clamp m>=1, avoids separate check
    else
        m = max(1, n - floor(theta));
    end

    %------------------
    % Decomposition
    %------------------
    % Precompute shared constants
    sqrtN_pi  = sqrt(N / pi);
    sqrtn_pi  = sqrt(n / pi);
    pi_over_N = pi / N;
    inv2m     = 1 / (2 * m);

    % STEP 1: DCT of every 3rd element (indices 2,5,8,...)
    ap     = a(2:3:end);                        % length n
    alphaf = dct(ap * pi_over_N) * sqrtn_pi;

    % STEP 2: DCT of a with those same elements zeroed
    app        = a;
    app(2:3:end) = 0;
    y2 = dct(app * pi_over_N) * sqrtN_pi;

    % Build betaf — first copy y2(1:n), then overwrite tail vectorized
    betaf = y2(1:n);

    % Transition-band indices
    k_range = (n-m+1 : n-1)';          % column vector, length m-1
    if ~isempty(k_range)
        w1 = (m + n - k_range) * inv2m;
        w2 = (m - n + k_range) * inv2m;
        betaf(k_range + 1) = w1 .* y2(k_range + 1) - w2 .* y2(2*n - k_range + 1);
    end

    % STEP 3: mu weights and combined coefficient ab
    mu = ones(n, 1);
    if ~isempty(k_range)
        mu(k_range + 1) = (m^2 + (n - k_range).^2) / (2 * m^2);
    end

    ab = (alphaf + betaf) ./ mu;

    % Scaling coefficients
    anew = dct(ab, 'Type', 3) * sqrtn_pi * fsca;

    % STEP 4: build sparse-like vector c, inverse DCT, subtract
    c = zeros(N, 1);
    c(1:n-m+1) = ab(1:n-m+1);

    if ~isempty(k_range)
        c(k_range + 1)         =  ab(k_range + 1) .* w1;
        c(2*n - k_range + 1)   = -ab(k_range + 1) .* w2;
    end

    y5           = dct(c, 'Type', 3) * sqrtN_pi;
    y5(2:3:end)  = 0;

    % Wavelet coefficients
    bnew         = (app - y5);
    bnew(2:3:end) = [];
    bnew         = bnew * fwav;
end