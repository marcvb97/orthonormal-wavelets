function [anew, bnew, nplus] = FT1step_m(a, theta)
% FT1step_m  One-step forward wavelet transform for a single vector.
%
% Handles the edge case where m = floor(theta*n) = 0 (e.g. theta = 0.1,
% n = 9 after several decomposition levels).
%
% INPUT:
%   a      - column vector to transform (length N1)
%   theta  - localization parameter:
%              theta < 1  -->  m = floor(theta * n)
%              theta >= 1 -->  m = n - floor(theta)
%
% OUTPUT:
%   anew   - scaling coefficients (column vector of length n)
%   bnew   - wavelet coefficients (column vector of length 2*n)
%   nplus  - number of elements appended to a to reach length N = 3*n
%              (appended elements equal the last entry of a)

%% Pad a to length N = 3*n
N1 = length(a);
r1 = rem(N1, 3);
if r1 == 0
    N     = N1;  n = N/3;
    nplus = 0;
else
    nplus = 3 - r1;
    for k = 1:nplus
        a(N1+k) = a(N1);
    end
    N = N1 + nplus;  n = N/3;
end

%% Compute localization degree m
if theta < 1
    m = floor(theta * n);
else
    m = n - floor(theta);
end
if m == 0           % guard against critical case (e.g. n=9, theta=0.1)
    m = 2;
end

%% Build weight vectors nu, nu3, mu, v
nu = ones(n, 1);
for r = n-m+1 : n-1
    nu(r+1) = (m^2 + (n-r)^2) / (2*m^2);
end

nu3 = ones(3*n, 1);
for r = 3*n-m+1 : 3*n-1
    nu3(r+1) = (m^2 + (3*n-r)^2) / (2*m^2);
end

mu = ones(n+m, 1);
for r = n-m+1 : n+m-1
    mu(r+1) = (m + n - r) / (2*m);
end

v      = zeros(3*n, 1);
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

%% Forward transform: upper (scaling) part
al  = dct(a, Type=2);
al  = al .* sqrt(1 ./ nu3);

al2 = zeros(n, 1);
for s = 0 : n-m
    al2(s+1) = al(s+1) / sqrt(nu(s+1));
end
for s = n-m+1 : n-1
    al2(s+1) = (al(s+1)*mu(s+1) - al(2*n-s+1)*mu(2*n-s+1)) / sqrt(nu(s+1));
end
ap = dct(al2, Type=3);

%% Forward transform: lower (wavelet) part
al3      = zeros(3*n, 1);
al3(0+1) = sqrt(2) / sqrt(3*v(2*n+1)) * al(2*n+1);
for s = 1 : n-m
    al3(s+1) = al(2*n-s+1) / sqrt(2*v(2*n-s+1));
end
for s = n-m+1 : n-1
    al3(s+1) = (mu(2*n-s+1)*al(s+1) + mu(s+1)*al(2*n-s+1)) / sqrt(2*v(2*n-s+1));
end
al3(n+1) = al(n+1);
for r = n+1 : n+m-1
    al3(r+1) = (mu(r+1)*al(2*n-r+1) + mu(2*n-r+1)*al(r+1)) / sqrt(2*v(r+1));
end
for r = n+m : 2*n-1
    al3(r+1) = al(r+1) / sqrt(2*v(r+1));
end
al3(2*n+1) = al(2*n+1) / sqrt(3*v(2*n+1));
for r = 2*n+1 : 3*n-m
    al3(r+1) = sqrt(3/2) * al(r+1) / sqrt(v(r+1));
end
for r = 3*n-m+1 : 3*n-1
    al3(r+1) = sqrt(3/2) * al(r+1) * sqrt(v(r+1));
end

bp          = dct(al3, 'Type', 3);
bp(2:3:end) = [];

anew = ap;
bnew = bp;

end
