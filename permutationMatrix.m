function [M] = permutationMatrix(vec)
  n = length(vec);
  M = zeros(n);
  for i = 1:n
    M(vec(i),i) = 1;
  end
end