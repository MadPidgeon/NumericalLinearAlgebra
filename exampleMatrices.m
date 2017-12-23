n = 20;
m = 10;



pseudoRandVec = ones(1,n); %numbers in [-504,504]
newNumber = 1;
generator = 46; %This is a generator in (Z/1019Z)*
somePrime = 1019; %This is a good prime: phi(1018) = 508.
for k = 1:n
  newNumber = mod(generator * newNumber, 1019);
  pseudoRandVec(k) = newNumber - 504;
end

%-------------diagonal matrices
D1 = diag(1:n); %positive definite matrix 
D2 = diag([ones(1,m), (m+1):n]);
D3 = diag(pseudoRandVec); %pseudorandom diagonal matrix, indefinite
I = diag([ones(1,n)]); %identity matrix

%-------------Bidiagonal matrices-------------------------------
a = 23; b = -29;

%with subdiagonal:
B1 = diag(a*ones(1,n),0) + diag(b*ones(1,n-1), -1);

%with superdiagonal:
B2 = diag(a*ones(1,n),0) + diag(b*ones(1,n-1), 1);

%------------Tridiagonal matrices--------------------------------
a = 8; b = 1; c = 9;
T1 = diag(a*ones(1,n), 0) + diag(b*ones(1,n-1), -1) + diag(c*ones(1,n-1), 1);

%------------Banded matrices------------------------------------------
lowerBorder = -5; upperBorder = 3;
bandWidth = upperBorder - lowerBorder + 1;
disp(bandWidth);
antiDiagBand = 1:bandWidth; %or adjust this to:
%antiDiagBand = antiDiagBand * antiDiagBand - 3 * antiDiagBand; %or someth. else

BD = zeros(n);
for k = lowerBorder:upperBorder
  diagLength = n - abs(k);
  factor = antiDiagBand(k-lowerBorder+1);
  BD =  BD + diag(factor * ones(1,diagLength), k);
  disp(k+1+lowerBorder);
end
%now BD (a BanDed matrix) is defined.

%---------------non-invertible matrix with one non-zero row--------------------
%you can add such a matrix to a banded matrix for example
i = 7; %7-th row non-zero
unitVec = zeros(n,1);
unitVec(i) = 1;
colVec = 1:n; %for example
oneNonZeroRow = unitVec * colVec;

%------------Haenkel matrices------------------------------------
%easy: flip a BD matrix (with flip(BD) for example).

%----------------pseudo-random matrix---------------------------
pseudoRandVec = ones(1,n*n); %numbers in [-500501,500501]
newNumber = 1;
generator = 46; %This is a generator in (Z/1001003Z)*
somePrime = 1001003; %This is a good prime: phi(1018) = 500500.
%500501 is a Sophie-Germain prime

for k = 1:(n*n)
  newNumber = mod(generator * newNumber, somePrime);
  pseudoRandVec(k) = newNumber - 500501;
end