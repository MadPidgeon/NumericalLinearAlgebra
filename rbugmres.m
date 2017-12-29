function [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = rbugmres( A, b, tol, true_x, x0, iter_count )
  if nargin < 5
    x0 = zeros( length( b ), 1 );
  end
  N = length(b);
  if nargin < 6
    iter_count = N;
  end
  n = 1;
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
  U_condition_numbers(n) = 1; % temp
  Acn = cond(A);
  Anrm = norm(A);
  Z = zeros(N);
  V = zeros(N);
  U = zeros(N);
  P = zeros(N);
  u = zeros(N,1);
  x = x0;
  alpha = zeros(N,1);
  beta = zeros(N,1);
  while residual_norms(n) > tol && n <= iter_count
    Z(:,n) = r / residual_norms(n);
    V(:,n) = A*Z(:,n);   
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
    % -----
    % viezigheid
    % -----
    backward_error(n) = norm(b-A*x)/(norm(x)*Anrm);
    forward_error(n) = norm(true_x-x)/norm(true_x);
    true_residual(n) = norm(b-A*x)/norm(b);
    updated_residual(n) = norm(r)/norm(b);
    Z_condition_numbers(n) = cond(Z(:,1:(n-1))); % mogelijk raar
    %Z_condition_numbers(n) = cond(Z(:,1:(n-1))'*Z(:,1:(n-1)));
    U_condition_numbers(n) = cond(U(:,1:(n-1)));
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
end

