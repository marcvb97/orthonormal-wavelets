function [anew, bnew, nplus] = FT1step_m_fast(a, theta)
% FT1step_m  One-step forward wavelet transform (optimized).

%% Pad a to length N = 3*n
N1 = length(a);
r1 = rem(N1, 3);
if r1 ~= 0
    nplus    = 3 - r1;
    a(end+1:end+nplus) = a(N1);   % vectorized padding
else
    nplus = 0;
end
N = length(a);  n = N / 3;

%% Compute localization degree m
if theta < 1
    m = max(2, floor(theta * n));  % clamp m>=2, absorbs m==0 guard
else
    m = max(2, n - floor(theta));
end

%% Precompute index ranges (reused across multiple sections)
inv2m  = 1 / (2*m);
inv2m2 = 1 / (2*m^2);

r_nu   = (n-m+1 : n-1)';        % length m-1
r_nu3  = (3*n-m+1 : 3*n-1)';
r_mu   = (n-m+1 : n+m-1)';
nr_mu  = n - r_mu;               % (n - r) for mu

%% Build weight vectors nu, nu3, mu, v  (fully vectorized)
nu          = ones(n, 1);
nu(r_nu+1)  = (m^2 + (n - r_nu).^2) * inv2m2;

nu3             = ones(3*n, 1);
nu3(r_nu3+1)    = (m^2 + (3*n - r_nu3).^2) * inv2m2;

mu          = ones(n+m, 1);
mu(r_mu+1)  = (m + nr_mu) * inv2m;   % (m + n - r)/(2m)

% Build v in four vectorised segments
v             = ones(3*n, 1);
v(1:n)        = 0;                    % indices 0..n-1 stay zero except v(n+1)=1 handled by ones above... 
% re-do cleanly:
v             = zeros(3*n, 1);
v(n+1)        = 1.0;
r_v1          = (n+1 : n+m-1)';
v(r_v1+1)     = (m^2 + (n - r_v1).^2) * inv2m2;
r_v2          = (n+m : 3*n-m)';
v(r_v2+1)     = 1.0;
r_v3          = (3*n-m+1 : 3*n-1)';
v(r_v3+1)     = (m^2 + (3*n - r_v3).^2) * inv2m2;

%% Precompute reciprocal square-roots (avoid repeated sqrt+divide)
sqrt_inv_nu3  = 1 ./ sqrt(nu3);       % length 3n
sqrt_inv_nu   = 1 ./ sqrt(nu);        % length n

%% Forward transform: upper (scaling) part
al  = dct(a, 'Type', 2) .* sqrt_inv_nu3;

% Vectorized al2 assembly
s0  = (0 : n-m)';                     % indices where formula is al/sqrt(nu)
s1  = (n-m+1 : n-1)';                 % transition band

al2             = zeros(n, 1);
al2(s0+1)       = al(s0+1) .* sqrt_inv_nu(s0+1);
mu_s1           = mu(s1+1);
mu_sym          = mu(2*n - s1 + 1);
al2(s1+1)       = (al(s1+1) .* mu_s1 - al(2*n - s1 + 1) .* mu_sym) ...
                  .* sqrt_inv_nu(s1+1);

anew = dct(al2, 'Type', 3);

%% Forward transform: lower (wavelet) part
% Precompute sqrt factors for v-indexed denominators
sqrt_inv_2v = 1 ./ sqrt(2 * v + (v == 0));   % guard /0 where v=0 (unused slots)
sqrt_3_over_v = sqrt(3 ./ max(v, eps));        % for the sqrt(3/2)/sqrt(v) segment

al3 = zeros(3*n, 1);

% s = 0
al3(1)      = sqrt(2) * sqrt_inv_2v(2*n+1) / sqrt(3/2) * al(2*n+1);
% simpler: 1/sqrt(3*v(2n+1)) * sqrt(2) --> keep original formula factored:
al3(1)      = al(2*n+1) * sqrt(2 / (3 * v(2*n+1)));

% s = 1..n-m
s_a         = (1 : n-m)';
al3(s_a+1)  = al(2*n - s_a + 1) .* sqrt_inv_2v(2*n - s_a + 1);

% s = n-m+1..n-1  (transition)
mu_s1_lo    = mu(2*n - s1 + 1);
mu_s1_hi    = mu(s1 + 1);
al3(s1+1)   = (mu_s1_lo .* al(s1+1) + mu_s1_hi .* al(2*n - s1 + 1)) ...
              .* sqrt_inv_2v(2*n - s1 + 1);

% s = n
al3(n+1)    = al(n+1);

% r = n+1..n+m-1
r_a         = (n+1 : n+m-1)';
mu_r        = mu(r_a + 1);
mu_sym_r    = mu(2*n - r_a + 1);
al3(r_a+1)  = (mu_r .* al(2*n - r_a + 1) + mu_sym_r .* al(r_a+1)) ...
              .* sqrt_inv_2v(r_a + 1);

% r = n+m..2n-1
r_b         = (n+m : 2*n-1)';
al3(r_b+1)  = al(r_b+1) .* sqrt_inv_2v(r_b+1);

% r = 2n
al3(2*n+1)  = al(2*n+1) / sqrt(3 * v(2*n+1));

% r = 2n+1..3n-m
r_c         = (2*n+1 : 3*n-m)';
al3(r_c+1)  = sqrt(3/2) * al(r_c+1) .* sqrt_inv_2v(r_c+1) * sqrt(2);
% simplify: sqrt(3/2)/sqrt(v) = sqrt(3/(2v)) --> use sqrt_3_over_v/sqrt(2)
al3(r_c+1)  = al(r_c+1) .* sqrt(3 ./ (2 * v(r_c+1)));

% r = 3n-m+1..3n-1
al3(r_v3+1) = sqrt(3/2) * al(r_v3+1) .* sqrt(v(r_v3+1));

%% Wavelet coefficients
bp          = dct(al3, 'Type', 3);
bp(2:3:end) = [];

anew = anew;   % already set above
bnew = bp;
end