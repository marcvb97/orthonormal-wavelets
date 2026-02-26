%% test FWTmatrix1 and IFWTmatrix
theta = 0.5;
kmax = 4;
dmin = 1;
A = rand(1000,1000);
[Adec, lrow, lcol] = FWTmatrix1_m(A, theta, kmax, dmin);
Arec = IFWTmatrix_m(Adec, lrow, lcol, theta);
norm(Arec(1:1000,1:1000)-A)