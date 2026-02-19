function A = TR2D(F,theta);
[nr,nc] = size(F);

nuc = compute_nu(nc,round(theta*nc)).';
for i = 1:nr
    F(i,:) = TR1D(F(i,:),nuc);
end

nur = compute_nu(nr,round(theta*nr));
for j = 1:nc
    F(:,j) = TR1D(F(:,j),nur);
end
A = F;
