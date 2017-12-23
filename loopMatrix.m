function [matrixList] = loopMatrix(n)
  matrixList = zeros(n,n,6);
  matrixList(:,:,1) = diag([ones(1,n)]);
  matrixList(:,:,2) = diag(1:n);
  m = floor(sqrt(n));
  matrixList(:,:,3) = diag([ones(1,m), (m+1):n]);
  
  %------------------generating randomness------------------------------------
  pseudoRandVec = ones(1,n*n); %numbers in [-500501,500501]
  newNumber = 1;
  generator = 46; %This is a generator in (Z/1001003Z)*
  somePrime = 1001003; %This is a good prime: phi(1018) = 500500.
  %500501 is a Sophie-Germain prime
  for k = 1:(n*n)
    newNumber = mod(generator * newNumber, somePrime);
    pseudoRandVec(k) = newNumber - 500501;
  end
  %--------------------done!--------------------------------------------------
  
  matrixList(:,:,4) = diag(pseudoRandVec(1:n));
  
  %bidiagonal:
  matrixList(:,:,5) = diag(ones(1,n),0) - diag(ones(1,n-1), -1);
  matrixList(:,:,6) = diag(1000*ones(1,n),0) + diag(ones(1,n-1),1);
end