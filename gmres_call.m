function [x,res] = gmres_call( A, b, tol )
  [x,flag,relres,iter,res] = gmres(A,b,length(b),tol);
end 