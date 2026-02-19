% test the functions TR2D and ITR2D

F = rand(1000,900);

theta = 0.7;
A = TR2D(F,theta);
F2 = ITR2D(A,theta);
norm(F-F2) / norm(F)