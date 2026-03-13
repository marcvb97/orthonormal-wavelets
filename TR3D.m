function A = TR3D(F, theta)
[nr, nc, ns] = size(F);

% Transform along dim 2 (rows): pass each row as column, store back as row
nuc = compute_nu(nc, round(theta * nc));
for k = 1:ns
    for i = 1:nr
        F(i, :, k) = TR1D(F(i, :, k).', nuc).';  % col in, col out, store as row
    end
end

% Transform along dim 1 (columns): already column vectors
nur = compute_nu(nr, round(theta * nr));
for k = 1:ns
    for j = 1:nc
        F(:, j, k) = TR1D(F(:, j, k), nur);       % col in, col out
    end
end

% Transform along dim 3 (slices): squeeze to column
nus = compute_nu(ns, round(theta * ns));
for i = 1:nr
    for j = 1:nc
        slice = squeeze(F(i, j, :));               % ns x 1 column
        F(i, j, :) = TR1D(slice, nus);             % col in, col out
    end
end

A = F;