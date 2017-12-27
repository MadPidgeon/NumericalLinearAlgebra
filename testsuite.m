function [] = testsuite()
  sizes = [128 ... #,1000,10000,100000,1000000
    ];
  funcs = { @gmres_call, @rbugmres, @sgmres, @rbsgmres_fast  };
  for si = 1:length(sizes)
    n = sizes(si);
    matrices = loopMatrix(n);
    #disp(matrices);
    b = ones(n,1);
    for mi = 1:length(matrices)
      hold off;
      #disp(funcs{2});
      subplot(4,4,mi);
      for fi = 1:length(funcs)
        func = funcs{fi};
        [x,res] = func( matrices{1,mi}, b, 1e-8 );
        #disp(res);
        semilogy( res );
        %plot(res);
        hold on;
      end
    end
  end
end