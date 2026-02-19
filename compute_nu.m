function nu = compute_nu(n,m)
for r = 0:n-m
    nu(r+1) = 1.0;
end
for r = n-m+1:n-1
    nu(r+1) = (m^2+(n-r)^2) / (2*m^2);
end
nu = nu(:);