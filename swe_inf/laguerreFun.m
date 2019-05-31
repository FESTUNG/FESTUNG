
function ret = laguerreFun(i, X, beta)
ret = exp(-(beta/2).*X).*laguerrePol(i,X,beta);
end % function
