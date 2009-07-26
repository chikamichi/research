disp 'test de préconditionneur'
diag1 = ones(5,1).*1e-10;
diag3 = ones(5,1).*2.4.*1e-17;
diag2 = ones(5,1).*5.*1e8;
K = spdiags([diag1 diag2 diag3], -1:1:1, 5, 5);
b = rand(5,1);

x1 = K\b

M = spdiags(repmat([10^mantexpnt(diag1) 10^mantexpnt(diag2) 10^mantexpnt(diag3)],5,1), -1:1:1, 5, 5);
%inv(M')*K
x2 = pcg(K,b,[],[],M)


