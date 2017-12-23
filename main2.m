n=10;
m=3;
A = diag([ones(1,m), (m+1):n]);
b = ones(n,1);
[x,residual_norms] = rbugmres(A,b,1e-8);
x
residual_norms
disp(norm(A*x-b));