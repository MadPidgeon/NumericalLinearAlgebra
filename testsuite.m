function [] = testsuite( func, b )
  sizes = [10,100,1000,10000,100000,1000000];
  funcs = [ rbugmres, sgmres, rbsgmres ];
  for si = 1:length(sizes)
    n = sizes[si];
    matrices = loopmatrix(n);
    for mi = 1:size(matrices)(3)
      hold off
      for fi = 1:size(funcs)
        func = funcs(fi);
        [x,res] = func( A, b, 1e-8 );
        semilogy( res );
        hold on
      end
    end
  end
end 