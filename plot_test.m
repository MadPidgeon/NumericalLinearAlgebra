function [] = plot_test()
  sizes = [128 ... #,1000,10000,100000,1000000
    ];
  funcs = { @rbsgmres };
  for si = 1:length(sizes)
    n = sizes(si);
    matrices = loopMatrix(n);
    #disp(matrices);
    true_x = ones(n,1);
    for mi = 1:length(matrices)
      b = matrices{1,mi}*true_x;
      hold off;
      #disp(funcs{2});
      subplot(4,4,mi);
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
        %plot(res);
      end
    end
  end
  