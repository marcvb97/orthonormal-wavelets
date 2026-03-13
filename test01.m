% Quick diagnostic - run this before fixing
F_test = rand(100, 100, 100);
theta = 0.5;
nr = 100; nc = 100; ns = 100;

nur = compute_nu(nr, round(theta * nr));
col = F_test(:, 1, 1);           % 100x1 column
fprintf('Input col size:       %dx%d\n', size(col, 1), size(col, 2));
result = TR1D(col.', nur);
fprintf('TR1D output size:     %dx%d\n', size(result, 1), size(result, 2));

row = F_test(1, :, 1);           % 1x100 row
result2 = TR1D(row, nur);
fprintf('TR1D row output size: %dx%d\n', size(result2, 1), size(result2, 2));