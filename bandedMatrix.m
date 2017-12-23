function [BD] = bandedMatrix(n,antiDiagBand,lowerBorder)
  bandWidth = length(antiDiagBand);
  upperBorder = lowerBorder + bandWidth - 1;

  BD = zeros(n);
  for k = lowerBorder:upperBorder
    diagLength = n - abs(k);
    factor = antiDiagBand(k-lowerBorder+1);
    BD =  BD + diag(factor * ones(1,diagLength), k);
    disp(k+1+lowerBorder);
  end
  %now BD (a BanDed matrix) is defined.
end