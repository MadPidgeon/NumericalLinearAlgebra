function [x,residual_norms] = rbsgmres_fast( A, b, tol, x0 )
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
  u = zeros(N,1);
  alpha = zeros(N,1);
  beta = zeros(N,1);
  while residual_norms(n) > tol && n <= N
    Z(:,n) = r / residual_norms(n);
    V(:,n) = A*Z(:,n);
    %for i = 1:(n-1)
    %  U(i,n) = dot(V(:,n),V(:,i));
    %  V(:,n) -= U(i,n)*V(:,i);
    %end
    %gamma = norm(V(:,n));
    %U(n,n) = gamma;
    %V(:,n) /= gamma;
    [V(:,n),U(1:n,n)] = mgorth(V(:,n), V(:,1:(n-1)));
    alpha(n) = dot(r,V(:,n));
    r -= alpha(n)*V(:,n);
    n += 1;
    residual_norms(n) = norm(r);
  end
  residual_norms = residual_norms(1:n);
  n -= 1;
  t = U(:,1:n) \ alpha;
  x = x0 + Z(:,1:n)*t;
end

