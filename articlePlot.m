function [] = articlePlot()
  [M0,M1,M2] = articleMatrices();
  index = 1;
  if( index == 1 )
    n = length(M1);
    true_x = ones(n,1);
    b = M1*true_x;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = rbsgmres( M1, b, 1e-8, true_x );
    semilogy( backward_error );
    hold on;
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    hold off;
  end
end