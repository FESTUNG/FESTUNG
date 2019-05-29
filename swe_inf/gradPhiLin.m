
function ret = gradPhiLin(i, m, X1, X2)
switch m
  case 1
    switch i
      case 1,  ret = zeros(size(X2));
      case 2,  ret = zeros(size(X1));
    end
  case 2
    switch i
      case 1,  ret = zeros(size(X2));       
      case 2,  ret = ones(size(X2));
    end % switch
end % switch
end % function
