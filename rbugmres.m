function [x,residual_norms] = rbugmres( A, b, tol, x0 )
  if nargin < 4
    x0 = zeros( length( b ), 1 );
  end
  N = length(b);
  n = 1;
  residual_norms = zeros(N,1);
  r = b - A*x0;
  residual_norms(n) = norm(r);
  Z = zeros(N);
  V = zeros(N);
  U = zeros(N);
  P = zeros(N);
  u = zeros(N,1);
  x = x0;
  alpha = zeros(N,1);
  beta = zeros(N,1);
  while residual_norms(n) > tol && n <= N
    Z(:,n) = r / residual_norms(n);
    V(:,n) = A*Z(:,n);   
    #for i = 1:(n-1)
    #  U(i,n) = dot(V(:,n),V(:,i));
    #  V(:,n) -= U(i,n)*V(:,i);
    #end
    #gamma = norm(V(:,n));
    #U(n,n) = gamma;
    #V(:,n) /= gamma;
    [V(:,n),U(1:n,n)] = mgorth(V(:,n), V(:,1:(n-1)));
    for j = 1:N
      s = P(j,1:(n-1))*U(1:(n-1),n);
      P(j,n) = ( Z(j,n) - s ) / U(n,n);
    end
    alpha(n) = dot(r,V(:,n));
    r -= alpha(n)*V(:,n);
    x += alpha(n)*P(:,n);
    n += 1;
    residual_norms(n) = norm(r);
  end
  residual_norms = residual_norms(1:n);
end

