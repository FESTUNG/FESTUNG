function ret = legendre01(i,X)
switch i
  case 1,  ret = ones(size(X));
  case 2,  ret = sqrt(12)   .*(                                X - 0.5 );
  case 3,  ret = sqrt(180)  .*(                     X.^2 -     X + 1/6 );
  case 4,  ret = sqrt(2800) .*(          X.^3 - 3/2*X.^2 + 3/5*X - 1/20);  
  case 5,  ret = sqrt(44100).*(X.^4 - 2.*X.^3 + 9/7*X.^2 - 2/7*X + 1/70);
end % switch
end % function