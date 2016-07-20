function ret = legendre01Prime(i,X)
switch i
  case 1,  ret =             zeros(size(X));
  case 2,  ret = sqrt(12)   .*ones(size(X));
  case 3,  ret = sqrt(180)  .*(                      2*X -   1);
  case 4,  ret = sqrt(2800) .*(          3*X.^2 -    3*X + 3/5);  
  case 5,  ret = sqrt(44100).*(4*X.^3 - 6.*X.^2 + 18/7*X - 2/7);
end % switch
end % function