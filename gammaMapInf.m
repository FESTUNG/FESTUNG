function [X1, X2] = gammaMapInf(n, S)
S = S(:)';
switch n
  case 1,  X1 = S;               X2 = ones(size(S));
  case 2,  X1 = zeros(size(S));  X2 = 1 - S; %endliche Kante!
  case 3,  X1 = S;               X2 = zeros(size(S));
end % switch
end % function
