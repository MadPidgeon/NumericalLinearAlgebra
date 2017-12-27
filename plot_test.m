function [] = plot_test()
  funcs = { @rbugmres };
  [M0,M1,M2] = articleMatrices();
  matrices = {M0,M1,M2};
  #disp(matrices);
  for mi = 1:length(matrices)
    n = length(matrices{mi});
    %true_x = vectors{1,mi};
    true_x = ones(n,1);
    b = matrices{mi}*true_x;
    hold off;
    #disp(funcs{2});
    subplot(1,3,mi);
    for fi = 1:length(funcs)
      func = funcs{fi};
      [x,residual_norms,backward_error,forward_error,true_residual,updated_residual,Z_condition_numbers] = func( matrices{1,mi}, b, 1e-8, true_x );
      #disp(res);
      semilogy( residual_norms );
      hold on;
      semilogy( backward_error );
      semilogy( forward_error );
      semilogy( true_residual );
      semilogy( updated_residual );
      semilogy( Z_condition_numbers );
      %disp( Z_condition_numbers );
      %plot(res);
    end
  end
end
  