function [x,residual_norms] = gsa( A, b, tol, x0 )
  if nargin < 4
    x0 = zeros( length( b ), 1 );
  end
  n = 1;
  while norm(r) > tol
    
    
    n += 1;
  end