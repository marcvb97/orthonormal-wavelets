function F = ITR2D(A,theta);
[nr,nc] = size(A);

F = A;
nur = compute_nu(nr,round(theta*nr));
for j = 1:nc
    F(:,j) = ITR1D(F(:,j),nur);
end

nuc = compute_nu(nc,round(theta*nc)).';
for i = 1:nr
    F(i,:) = ITR1D(F(i,:),nuc);
end
