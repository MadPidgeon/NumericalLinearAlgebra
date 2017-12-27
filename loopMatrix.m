function [matrixList] = loopMatrix(n)
  matrixList = cell(1,16);
  matrixList(1) = sparse(diag([ones(1,n)]));
  matrixList(2) = sparse(diag(1:n));
  m = floor(sqrt(n));
  matrixList(3) = sparse(diag([ones(1,m), (m+1):n]));
  
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
  
  matrixList(4) = sparse(diag(pseudoRandVec(1:n)));
  
  %bidiagonal:
  matrixList(5) = sparse(bandedMatrix(n,[-1,1],-1));
  matrixList(6) = sparse(bandedMatrix(n,[1000,1],0));
  matrixList(7) = sparse(bandedMatrix(n,[1,1000],-1));
  
  %tridiagonal
  matrixList(8) = sparse(bandedMatrix(n,[-1,2,-1],-1));
  matrixList(9) = sparse(bandedMatrix(n,[1,1000,-1],-1));
  matrixList(10) = sparse(bandedMatrix(n,[-1,1,3],-1));
  
  %banded matrix:
  matrixList(11) = sparse(bandedMatrix(n,1:m,0));
  matrixList(12) = sparse(bandedMatrix(n,1:m,-floor(m/2)));
  
  unitVec = zeros(n,1);
  unitVec(1) = 1;
  colVec = 1:n; %for example
  oneNonZeroRow = sparse(unitVec * colVec);
  
  %a combination of a tridiagonal matrix with a first-row-non-zero matrix:
  matrixList(13) = sparse(matrixList{8} + oneNonZeroRow);
  
  R = zeros(n); %this will be the random matrix
  for i = 1:n
    for j = 1:n
      R(i,j) = pseudoRandVec(n*(i-1)+j);
    end
  end
  matrixList(14) = R;
  
  matrixList(15) = hadamard(n - mod(n,4));
  
  temp = diag(1:n);
  temp(1) = 1.0e-008;
  temp(2) = 2.0e-008;
  theta = pi/4;
  temp = givensRotation(n,1,m,theta) * temp * givensRotation(n,1,n,-theta);
  matrixList(16) = sparse(temp);
  
end