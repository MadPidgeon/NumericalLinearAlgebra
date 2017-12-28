function [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = sgmres( A, b, tol, true_x, x0, iter_count )
  if nargin < 5
    x0 = zeros( length( b ), 1 );
  end
  if nargin < 6
    iter_count = N;
  end
  N = length(b);
  n = 1;
  Acn = cond(A);
  Anrm = norm(A);
  residual_norms = zeros(N,1);
  backward_error = zeros(N,1);
  forward_error = zeros(N,1);
  true_residual = zeros(N,1);
  updated_residual = zeros(N,1);
  Z_condition_numbers = zeros(N,1);
  U_condition_numbers = zeros(N,1);
  r = b - A*x0;
  residual_norms(n) = norm(r);
  backward_error(n) = 1; % temp
  forward_error(n) = 1; % temp
  true_residual(n) = norm(b-A*x0)/norm(b);
  updated_residual(n) = norm(r)/norm(b);
  Z_condition_numbers(n) = 1; % temp
  U_condition_numbers(n) = 1;
  Z = zeros(N);
  V = zeros(N);
  U = zeros(N);
  u = zeros(N,1);
  alpha = zeros(N,1);
  beta = zeros(N,1);
  Z(:,n) = r / residual_norms(n);
  while residual_norms(n) > tol && n <= iter_count
    V(:,n) = A*Z(:,n);
    [V(:,n),U(1:n,n)] = mgorth(V(:,n), V(:,1:(n-1)));
    alpha(n) = dot(r,V(:,n));
    r -= alpha(n)*V(:,n);
    n += 1;
    residual_norms(n) = norm(r);
    Z(:,n) = V(:,n-1);
    % -----
    % viezigheid
    % -----
    x = x0 + Z(:,1:n-1)*(U(:,1:n-1) \ alpha);
    backward_error(n) = norm(b-A*x)/(norm(x)*Anrm);
    forward_error(n) = norm(true_x-x)/norm(true_x);
    true_residual(n) = norm(b-A*x)/norm(b);
    updated_residual(n) = norm(r)/norm(b);
    Z_condition_numbers(n) = cond(Z(:,1:(n-1))); % mogelijk raar
    U_condition_numbers(n) = cond(U(1:(n-1),1:(n-1)));
    % -----
    % end  
    % -----
  end
  residual_norms = residual_norms(1:n);
  backward_error = backward_error(1:n);
  forward_error = forward_error(1:n);
  true_residual = true_residual(1:n);
  updated_residual = updated_residual(1:n);
  Z_condition_numbers = Z_condition_numbers(1:n);
  U_condition_numbers = U_condition_numbers(1:n);
  n -= 1;
  t = U(:,1:n) \ alpha;
  x = x0 + Z(:,1:n)*t;
end

