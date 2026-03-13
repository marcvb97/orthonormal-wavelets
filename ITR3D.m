function F = ITR3D(A, theta)
[nr, nc, ns] = size(A);
F = A;

% Inverse transform along dim 3 (slices)
nus = compute_nu(ns, round(theta * ns));
for i = 1:nr
    for j = 1:nc
        slice = squeeze(F(i, j, :));               % ns x 1 column
        F(i, j, :) = ITR1D(slice, nus);            % col in, col out
    end
end

% Inverse transform along dim 1 (columns)
nur = compute_nu(nr, round(theta * nr));
for k = 1:ns
    for j = 1:nc
        F(:, j, k) = ITR1D(F(:, j, k), nur);      % col in, col out
    end
end

% Inverse transform along dim 2 (rows)
nuc = compute_nu(nc, round(theta * nc));
for k = 1:ns
    for i = 1:nr
        F(i, :, k) = ITR1D(F(i, :, k).', nuc).';  % col in, col out, store as row
    end
end