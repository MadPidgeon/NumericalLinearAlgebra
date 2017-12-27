function [] = testsuite2()
  sizes = [100 ... #,1000,10000,100000,1000000
    ];
  funcs = { @rbugmres, @sgmres, @rbsgmres_fast, @gmres_call };
  for si = 1:length(sizes)
    n = sizes(si);
    [matrices, vectors] = loopMatrix(n); %matrices = A, vectors = x, in Ax=b.
    #disp(matrices);
    for mi = 1:length(matrices)
      disp("MATRIX:");
      disp(mi);
      hold off;
      #disp(funcs{2});
      subplot(6,4,mi);
      for fi = 1:length(funcs)
        disp(fi);
        fflush(stdout);
        func = funcs{fi};
        b = matrices{1,mi} * vectors{1,mi};
        [x,res] = func( matrices{1,mi}, b, 1e-8 );
        %disp(res);
        semilogy( res );
        %plot(res);
        hold on;
      end
    end
  end
end