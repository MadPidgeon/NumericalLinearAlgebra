function [] = articlePlot()
  [M0,M1,M2] = articleMatrices();
  index = 0;
  if( index == 0 )
    n = length(M0);
    b = ones(n,1);
    true_x = M0\b;
    iter_count = 100;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = generalized_simpler_approach( M0, b, 1e-8, true_x, zeros(n,1), iter_count );
    hold off;
    semilogy( backward_error );
    hold on;
    semilogy( forward_error );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = generalized_update_approach( M0, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    hold on;
    semilogy( forward_error );
    hold off;
    legend('Location','southwest');
    legend('boxoff');
    legend('backward error (simpler)', 'forward error (simpler)','backward error (updated)', 'forward error (updated)');
  end
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
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = sgmres_explict_mgs( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( updated_residual );
    semilogy( true_residual );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = orthodir( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    hold off;
    legend('Location','southwest');
    legend('boxoff');
    legend('Simpler GMRES backward error','Simpler GMRES true residual','Simpler GMRES residual approximation','Simpler GMRES forward error','bad SGMRES residual approximation','bad SGMRES true residual','ORTHODIR backward error','ORTHODIR true residual','ORTHODIR residual approximation','ORTHODIR forward error');
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
  if( index == 3 )
    n = length(M1);
    true_x = ones(n,1);
    b = M1*true_x;
    iter_count = 60;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = rbsgmres( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( eps*Z_condition_numbers );
    hold on;
    semilogy( eps*U_condition_numbers );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = sgmres( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( eps*Z_condition_numbers );
    semilogy( eps*U_condition_numbers );
    %[x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = sgmres_explict_mgs( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    %semilogy( eps*U_condition_numbers );
    %[x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = orthodir( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    %semilogy( eps*U_condition_numbers );
    base = eps*cond(M1)*ones(length(Z_condition_numbers),1);
    semilogy( base );
    hold off;
    legend('Location','northwest');
    legend('boxoff');
    legend('RBSGMRES condition number tildeR','RBSGMRES condition number U','Simpler GMRES condition number tildeR','Simpler GMRES condition number U','\epsilon\kappa(A)');
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
  if( index == 5 )
    n = length(M2);
    true_x = ones(n,1);
    b = M2*true_x;
    iter_count = 220;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = rbsgmres( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    hold on;
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = gcr( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( backward_error );
    semilogy( true_residual );
    semilogy( updated_residual );
    semilogy( forward_error );
    hold off;
    legend('Location','southwest');
    legend('boxoff');
    legend('RBSGMRES backward error','RBSGMRES true residual','RBSGMRES residual approximation','RBSGMRES forward error','GCR backward error','GCR true residual','GCR residual approximation','GCR forward error');
  end
  if( index == 6 )
    n = length(M2);
    true_x = ones(n,1);
    b = M2*true_x;
    iter_count = 220;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = rbsgmres( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( eps*Z_condition_numbers );
    hold on;
    semilogy( eps*U_condition_numbers );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = sgmres( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( eps*Z_condition_numbers );
    semilogy( eps*U_condition_numbers );
    base = eps*cond(M2)*ones(length(Z_condition_numbers),1);
    semilogy( base );
    hold off;
    legend('Location','northwest');
    legend('boxoff');
    legend('RBSGMRES condition number tildeR','RBSGMRES condition number U','Simpler GMRES condition number tildeR','Simpler GMRES condition number U','\epsilon\kappa(A)');
  end
end