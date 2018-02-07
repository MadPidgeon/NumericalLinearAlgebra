% The RB-GMRES algorithm modified to compute and return various statistics about the convergence assuming the exact solution is known
% INPUT:
%   A          The matrix in question
%   b          The image
%   tol        The tolerance in residual norm, i.e. the value such that if norm(r)<tol, the algorithm terminates
%   true_x     The value of x solving A*x=b exactly
%   x0         The initial guess of x to start with (is 0 by default)
%   iter_count An upper bound on the number of iteration to run (is length(b) by default)
% OUTPUT:
%   x                   The computed solution to A*x=b
%   residual_norms      The value of norm(r) in each iteration
%   backward_error      The value of norm(b-A*x)/(norm(x)*norm(A)) in each iteration
%   forward_error       The value of norm(true_x-x)/norm(true_x) in each iteration
%   true_residual       The value of norm(b-A*x)/norm(b) in each iteration
%   updated_residual    The value of norm(r)/norm(b) in each iteration
%   Z_condition_numbers The condition number of Z = [ r_n / norm(r_n) ]_n in each iteration
%   U_condition_numbers The condition number of U in each iteration

function [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = rbsgmres( A, b, tol, true_x, x0, iter_count )
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
  r = b - A*x0;
  residual_norms(n) = norm(r);
  backward_error(n) = 1;
  forward_error(n) = 1;
  true_residual(n) = norm(b-A*x0)/norm(b);
  updated_residual(n) = norm(r)/norm(b);
  Z_condition_numbers(n) = 1;
  U_condition_numbers(n) = 1;
  Acn = cond(A);
  Anrm = norm(A);
  Z = zeros(N);
  V = zeros(N);
  U = zeros(N);
  u = zeros(N,1);
  alpha = zeros(N,1);
  beta = zeros(N,1);
  while residual_norms(n) > tol && n <= iter_count
    Z(:,n) = r / residual_norms(n);
    V(:,n) = A*Z(:,n);
    [V(:,n),U(1:n,n)] = mgorth(V(:,n), V(:,1:(n-1)));
    alpha(n) = dot(r,V(:,n));
    r -= alpha(n)*V(:,n);
    n += 1;
    residual_norms(n) = norm(r);
    x = x0 + Z(:,1:n-1)*(U(:,1:n-1) \ alpha);
    backward_error(n) = norm(b-A*x)/(norm(x)*Anrm);
    forward_error(n) = norm(true_x-x)/norm(true_x);
    true_residual(n) = norm(b-A*x)/norm(b);
    updated_residual(n) = norm(r)/norm(b);
    Z_condition_numbers(n) = cond(Z(:,1:(n-1)));
    U_condition_numbers(n) = cond(U(:,1:(n-1)));
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

