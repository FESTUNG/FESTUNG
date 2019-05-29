function ret = laguerrePolPrime(i, X, beta)
ret = zeros(size(X));
for j = 1:(i-1)
  ret = ret + laguerrePol(j, X, beta);
end % for
ret = -beta*ret;
end % function