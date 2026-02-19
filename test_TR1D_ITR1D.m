% test the functions TR1D and ITR1D
n = 1000;
m = 100;
nu = compute_nu(n,m);
f = randn(n,1);
a = TR1D(f,nu);
f2 = ITR1D(a,nu);
norm(f-f2)/norm(f)