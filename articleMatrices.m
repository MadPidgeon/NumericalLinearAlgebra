function [A,B,C] = articleMatrices()
  %----------Figure 2.1 matrix----- = G_1 D G_2^T --------------------------
  temp = diag(1:100);
  temp(1) = 1.0e-008;
  temp(2) = 2.0e-008;
  theta = pi/4;
  temp = givensRotation(100,1,10,theta)*temp*givensRotation(100,1,100,theta)';
  A = sparse(temp); %matrix of figure 2.1
  
  %---------- FS 183 6 (from matrix market) -----------------------------------
  B = mmread('fs_183_6.mtx');
  
  %-------- STEAM1 (from matrix market) -------------------------------------
  C = mmread('steam1.mtx');
end