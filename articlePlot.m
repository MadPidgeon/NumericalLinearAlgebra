function [] = articlePlot()
  [M0,M1,M2] = articleMatrices();
  index = 3;
  if( index == 0 )
    n = length(M0);
    b = ones(n,1);
    true_x = M0\b;
    iter_count = 100;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = generalized_simpler_approach( M0, b, 1e-8, true_x, zeros(n,1), iter_count );
    hold off;
    semilogy( 0:iter_count, backward_error, "-", "linewidth", 2, "color", "k" );
    hold on;
    semilogy( 0:iter_count, forward_error, "--", "linewidth", 2, "color", "k" );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = generalized_update_approach( M0, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( 0:iter_count, backward_error,"-", "linewidth", 1, "color", "k" );
    hold on;
    semilogy( 0:iter_count, forward_error,"--", "linewidth", 1, "color", "k" );
    hold off;
    axis([0 iter_count]);
    legend('Location','southwest');
    legend('boxoff');
    legend('backward error (simpler)', 'forward error    (simpler)','backward error (updated)', 'forward error    (updated)');
  end
  if( index == 1 )
    n = length(M1);
    true_x = ones(n,1);
    b = M1*true_x;
    iter_count = 60;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = sgmres( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( 0:iter_count, backward_error, "-", "linewidth", 2, "color", "k" );
    hold on;
    semilogy( 0:iter_count, true_residual, "--", "linewidth", 2, "color", "k" );
    semilogy(  0:iter_count,forward_error, "-.", "linewidth", 2, "color", "k" );
    semilogy(  0:iter_count,updated_residual, ":", "linewidth", 1, "color", "k" );
    %[x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = sgmres_explict_mgs( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    %semilogy( updated_residual );
    %semilogy( true_residual );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = generalized_update_orthodir( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( 0:iter_count, backward_error, "-", "linewidth", 1, "color", "k" );
    semilogy( 0:iter_count, true_residual, "--", "linewidth", 1, "color", "k" );
    %semilogy( updated_residual );
    semilogy( 0:iter_count, forward_error, "-.", "linewidth", 1, "color", "k" );
    hold off;
    axis([0 iter_count]);
    legend('Location','southoutside');
    legend('boxoff');
    legend('Simpler GMRES backward error','Simpler GMRES relative residual norm','Simpler GMRES relative error norm','relative updated residual norm','ORTHODIR backward error','ORTHODIR relative residual norm','ORTHODIR relative error norm');
  end
  if( index == 2 )
    n = length(M1);
    true_x = ones(n,1);
    b = M1*true_x;
    iter_count = 60;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = rbsgmres( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    t = length(backward_error)-1;
    semilogy( 0:t, backward_error, "-", "linewidth", 2, "color", "k" );
    hold on;
    semilogy( 0:t, true_residual, "--", "linewidth", 2, "color", "k" );
    semilogy(  0:t,forward_error, "-.", "linewidth", 2, "color", "k" );
    semilogy(  0:t,updated_residual, ":", "linewidth", 1, "color", "k" );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = gcr( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    t = length(backward_error)-1;    
    semilogy( 0:t, backward_error, "-", "linewidth", 1, "color", "k" );
    semilogy( 0:t, true_residual, "--", "linewidth", 1, "color", "k" );
    semilogy( 0:t, forward_error, "-.", "linewidth", 1, "color", "k" );
    hold off;
    axis([0 t]);
    legend('Location','southoutside');
    legend('boxoff');
    legend('RBSGMRES backward error','RBSGMRES relative residual norm','RBSGMRES relative error norm','relative updated residual norm','GCR backward error','GCR relative residual norm','GCR relative error norm');
  end
  if( index == 3 )
    n = length(M1);
    true_x = ones(n,1);
    b = M1*true_x;
    iter_count = 60;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = rbsgmres( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    t = length(forward_error) - 1;
    semilogy( 0:t, eps*Z_condition_numbers, "-", "linewidth", 2, "color", "k" );
    hold on;
    semilogy( 0:t, eps*U_condition_numbers, "--", "linewidth", 2, "color", "k" );
    % disp(Z_condition_numbers); last is infty
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = sgmres( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    t = length(forward_error) - 1;
    semilogy( -1:(t-1), eps*Z_condition_numbers, "-", "linewidth", 1, "color", "k" );
    semilogy( -1:(t-1), eps*U_condition_numbers, "--", "linewidth", 1, "color", "k" );
    base = eps*cond(M1)*ones(length(Z_condition_numbers),1);
    semilogy( 0:t, base, "-.", "linewidth", 1, "color", "k" );
    hold off;
    axis([0 t]);
    legend('Location','southoutside');
    legend('boxoff');
    legend('RBSGMRES condition number R=\tilde{R}','RBSGMRES condition number U','Simpler GMRES condition number R=V','Simpler GMRES condition number U','\epsilon\kappa(A)');
  end
  if( index == 4 )
    n = length(M2);
    true_x = ones(n,1);
    b = M2*true_x;
    iter_count = 220;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = sgmres( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( 0:iter_count, backward_error, "-", "linewidth", 2, "color", "k" );
    hold on;
    semilogy( 0:iter_count, true_residual, "--", "linewidth", 2, "color", "k" );
    semilogy(  0:iter_count,forward_error, "-.", "linewidth", 2, "color", "k" );
    semilogy(  0:iter_count,updated_residual, ":", "linewidth", 1, "color", "k" );
    %[x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = sgmres_explict_mgs( M1, b, 1e-8, true_x, zeros(n,1), iter_count );
    %semilogy( updated_residual );
    %semilogy( true_residual );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = generalized_update_orthodir( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( 0:iter_count, backward_error, "-", "linewidth", 1, "color", "k" );
    semilogy( 0:iter_count, true_residual, "--", "linewidth", 1, "color", "k" );
    %semilogy( updated_residual );
    semilogy( 0:iter_count, forward_error, "-.", "linewidth", 1, "color", "k" );
    hold off;
    axis([0 iter_count]);
    legend('Location','southoutside');
    legend('boxoff');
    legend('Simpler GMRES backward error','Simpler GMRES relative residual norm','Simpler GMRES relative error norm','relative updated residual norm','ORTHODIR backward error','ORTHODIR relative residual norm','ORTHODIR relative error norm');
  end
  if( index == 5 )
    n = length(M2);
    true_x = ones(n,1);
    b = M2*true_x;
    iter_count = 220;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = rbsgmres( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    t = length(backward_error)-1;
    semilogy( 0:t, backward_error, "-", "linewidth", 2, "color", "k" );
    hold on;
    semilogy( 0:t, true_residual, "--", "linewidth", 2, "color", "k" );
    semilogy(  0:t,forward_error, "-.", "linewidth", 2, "color", "k" );
    semilogy(  0:t,updated_residual, ":", "linewidth", 1, "color", "k" );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = gcr( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    t = length(forward_error) - 1;
    semilogy( -1:(t-1), eps*Z_condition_numbers, "-", "linewidth", 1, "color", "k" );
    semilogy( -1:(t-1), eps*U_condition_numbers, "--", "linewidth", 1, "color", "k" );
    base = eps*cond(M1)*ones(length(Z_condition_numbers),1);
    semilogy( 0:t, base, "-.", "linewidth", 1, "color", "k" );
    hold off;
    axis([0 t]);
    legend('Location','southoutside');
    legend('boxoff');
    legend('RBSGMRES backward error','RBSGMRES relative residual norm','RBSGMRES relative error norm','relative updated residual norm','GCR backward error','GCR relative residual norm','GCR relative error norm');
  end
  if( index == 6 )
    n = length(M2);
    true_x = ones(n,1);
    b = M2*true_x;
    iter_count = 220;
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = rbsgmres( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    t = length(forward_error) - 1;
    semilogy( 0:t, eps*Z_condition_numbers, "-", "linewidth", 2, "color", "k" );
    hold on;
    semilogy( 0:t, eps*U_condition_numbers, "--", "linewidth", 2, "color", "k" );
    [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers,U_condition_numbers] = sgmres( M2, b, 1e-8, true_x, zeros(n,1), iter_count );
    semilogy( -1:(t-1), eps*Z_condition_numbers, "-", "linewidth", 1, "color", "k" );
    semilogy( -1:(t-1), eps*U_condition_numbers, "--", "linewidth", 1, "color", "k" );
    base = eps*cond(M1)*ones(length(Z_condition_numbers),1);
    semilogy( 0:t, base, "-.", "linewidth", 1, "color", "k" );
    axis([0 t]);
    legend('Location','southoutside');
    legend('boxoff');
    legend('RBSGMRES condition number R','RBSGMRES condition number U','Simpler GMRES condition number R','Simpler GMRES condition number U','\epsilon\kappa(A)');
  end
end