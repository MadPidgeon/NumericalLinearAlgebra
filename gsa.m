function [x,residual_norms] = gsa( A, b, tol, x0 )
  if nargin < 4
    x0 = zeros( length( b ), 1 );
  end
  N = length(b);
  n = 1;
  residual_norms = zeros(N,1);
  r = b - A*x0;
  Z = zeros(N);
  while norm(r) > tol
    
    
    
    residual_norms[n] = norm(r);
    n += 1;
  end