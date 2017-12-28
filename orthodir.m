function [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = orthodir( A, b, tol, true_x, x0, iter_count )
  if nargin < 5
    x0 = zeros( length( b ), 1 );
  end
  if nargin < 6
    iter_count = N;
  end
  N = length(b);
  n = 1;
  residual_norms = zeros(N,1);
  backward_error = zeros(N,1);
  forward_error = zeros(N,1);
  true_residual = zeros(N,1);
  updated_residual = zeros(N,1);
  Z_condition_numbers = ones(N,1);
  Acn = cond(A);
  Anrm = norm(A);
  r = b - A*x0;
  residual_norms(n) = norm(r);
  backward_error(n) = 1; % temp
  forward_error(n) = 1; % temp
  true_residual(n) = norm(b-A*x0)/norm(b);
  updated_residual(n) = norm(r)/norm(b);
  Z_condition_numbers(n) = 1; % temp
  x = x0;
  U = zeros(N);
  C = zeros(N);
  sigma = zeros(N,1);
  alpha = zeros(N,1);
  U(:,n) = r;
  while residual_norms(n) > tol && n <= iter_count
    U(:,n) = r;
    C(:,n) = A*U(:,n);
    for j = 1:(n-1)
      beta = dot(C(:,j),C(:,n)) / sigma(j);
      U(:,n) -= beta * U(:,j);
      C(:,n) -= beta * C(:,j);
    end
    sigma(n) = dot(C(:,n),C(:,n));
    alpha(n) = dot(C(:,n),r) / sigma(n);
    x += alpha(n)*U(:,n);
    r -= alpha(n)*C(:,n);
    n += 1;
    U(:,n) = C(:,n-1);
    residual_norms(n) = norm(r);
    % lelijke stuff
    backward_error(n) = norm(b-A*x)/(norm(x)*Anrm);
    forward_error(n) = norm(true_x-x)/norm(true_x);
    true_residual(n) = norm(b-A*x)/norm(b);
    updated_residual(n) = norm(r)/norm(b);
  end
  residual_norms = residual_norms(1:n);
  backward_error = backward_error(1:n);
  forward_error = forward_error(1:n);
  true_residual = true_residual(1:n);
  updated_residual = updated_residual(1:n);
  Z_condition_numbers = Z_condition_numbers(1:n);
end