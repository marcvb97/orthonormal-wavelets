function [aa, nminus] = IFT1step_m(a, b, theta)
% IFT1step_m  One-step inverse wavelet transform for a single vector.
%
% Handles the edge case where m = floor(theta*n) = 0 (e.g. theta = 0.1,
% n = 9 after several decomposition levels).
% If length(a) > n, the extra entries are discarded so that
% length(aa) = 3*n.
%
% INPUT:
%   a      - scaling coefficients (column vector)
%   b      - wavelet coefficients (column vector of even length 2*n)
%   theta  - localization parameter:
%              theta < 1  -->  m = floor(theta * n)
%              theta >= 1 -->  m = n - floor(theta)
%
% OUTPUT:
%   aa     - reconstructed column vector of length 3*n
%   nminus - number of entries discarded from a ( = length(a) - n )

%% Sizes and localization degree
n      = length(b) / 2;
nminus = length(a) - n;
N      = 3 * n;

if theta < 1
    m = floor(theta * n);
else
    m = n - floor(theta);
end
if m == 0           % guard against critical case (e.g. n=9, theta=0.1)
    m = 2;
end

%% Truncate scaling coefficients to length n
a = a(1:n);
b = b(1:2*n);

%% Step 1: alpha coefficients
al = dct(a, Type=2);

nu = ones(n, 1);
for r = n-m+1 : n-1
    nu(r+1) = (m^2 + (n-r)^2) / (2*m^2);
end
al = al ./ sqrt(nu);

%% Step 2: beta coefficients (interleave wavelet vector, then DCT)
u          = zeros(N, 1);
u(1:3:end) = b(1:2:end);
u(3:3:end) = b(2:2:end);
w          = dct(u, Type=2);

bet      = zeros(N, 1);
bet(n+1) = w(n+1);
for r = n+1 : 2*n-1
    bet(r+1) = (w(r+1) + w(2*n-r+1)) / sqrt(2);
end
bet(2*n+1) = (w(2*n+1) + sqrt(2)*w(1)) / sqrt(3);
for r = 2*n+1 : 3*n-1
    bet(r+1) = sqrt(3/2) * w(r+1);
end

v      = zeros(N, 1);
v(n+1) = 1.0;
for r = n+1 : n+m-1
    v(r+1) = (m^2 + (n-r)^2) / (2*m^2);
end
for r = n+m : 3*n-m
    v(r+1) = 1.0;
end
for r = 3*n-m+1 : 3*n-1
    v(r+1) = (m^2 + (3*n-r)^2) / (2*m^2);
end

bet = bet ./ sqrt(v);

%% Step 3: combine alpha and beta into gamma coefficients
mu = ones(n+m, 1);
for r = n-m+1 : n+m-1
    mu(r+1) = (m + n - r) / (2*m);
end

gam = zeros(N, 1);
for s = 0 : n-m
    gam(s+1) = al(s+1);
end
for s = n-m+1 : n-1
    gam(s+1) = al(s+1)*mu(s+1) + mu(2*n-s+1)*bet(2*n-s+1);
end
gam(n+1) = bet(n+1);
for r = n+1 : n+m-1
    gam(r+1) = bet(r+1)*mu(2*n-r+1) - mu(r+1)*al(2*n-r+1);
end
for r = n+m : 3*n-m
    gam(r+1) = bet(r+1);
end
for r = 3*n-m+1 : 3*n-1
    gam(r+1) = bet(r+1) * v(r+1);
end

%% Step 4: apply nu3 weighting and inverse DCT
nu3 = ones(N, 1);
for r = 3*n-m+1 : 3*n-1
    nu3(r+1) = (m^2 + (3*n-r)^2) / (2*m^2);
end

gam = gam ./ sqrt(nu3);
aa  = dct(gam, Type=3);

end
