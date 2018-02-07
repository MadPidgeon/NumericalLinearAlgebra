% Compute the Givens rotation as an nxn matrix rotating the i-th and j-th axes with angle theta 

function [M] = givensRotation(n, i, j, theta)
  M = diag(ones(1,n));
  M(i,i) = cos(theta);
  M(j,j) = cos(theta);
  M(i,j) = -sin(theta);
  M(j,i) = sin(theta);
  M = sparse(M);
end
