function [] = articlePlot()
  [M0,M1,M2] = articleMatrices();
  index = 4;
  if( index == 1 )
    n = length(M1);
    true_x = ones(n,1);
    b = M1*true_x;
    iter_count = 60;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = sgmres( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    hold on;
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = orthodir( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    hold off;
    legend('Location','southwest');
    legend('boxoff');
    legend('Simpler GMRES backward error','Simpler GMRES true residual','Simpler GMRES residual approximation','Simpler GMRES forward error','ORTHODIR backward error','ORTHODIR true residual','ORTHODIR residual approximation','ORTHODIR forward error');
  end
  if( index == 2 )
    n = length(M1);
    true_x = ones(n,1);
    b = M1*true_x;
    iter_count = 60;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = rbsgmres( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    hold on;
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = gcr( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    hold off;
    legend('Location','southwest');
    legend('boxoff');
    legend('RBSGMRES backward error','RBSGMRES true residual','RBSGMRES residual approximation','RBSGMRES forward error','GCR backward error','GCR true residual','GCR residual approximation','GCR forward error');
  end
  if( index == 4 )
    n = length(M2);
    true_x = ones(n,1);
    b = M2*true_x;
    iter_count = 220;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = sgmres( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    hold on;
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = orthodir( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    hold off;
    legend('Location','southwest');
    legend('boxoff');
    legend('Simpler GMRES backward error','Simpler GMRES true residual','Simpler GMRES residual approximation','Simpler GMRES forward error','ORTHODIR backward error','ORTHODIR true residual','ORTHODIR residual approximation','ORTHODIR forward error');
  end
end