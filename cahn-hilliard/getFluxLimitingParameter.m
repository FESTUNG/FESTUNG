function alphaEface = getFluxLimitingParameter(sysY, gamma, markE0TE0T, globFe, K, N, areaT)

Cmax = +1; Cmin = -1; epsilon = 1e-8;

cMean = sqrt(2) * sysY(1:N:K*N);
phiVec = sysY(K*N+1:end);

intFe = zeros(K,3);
alphaEface = zeros(K,3);

for n = 1 : 3
  intFe(:,n) = (globFe{n} * phiVec) ./ (sqrt(2) * areaT);
end % for n

alphaEplus = min( 1, gamma .* max(0, Cmax - cMean) ./ ( -sum( intFe .* (intFe < 0) , 2 ) + epsilon) );
alphaEminus = min( 1, gamma .* max(0, cMean - Cmin) ./ ( sum( intFe .* (intFe > 0) , 2 ) + epsilon) );

assert(all((0 <= alphaEplus) & (alphaEplus <= 1)), 'alphaEplus: %.3e %.3e', [max(max(alphaEplus)) - 1, min(min(alphaEplus))])
assert(all((0 <= alphaEminus) & (alphaEminus <= 1)), 'alphaEminus: %.3e %.3e', [max(max(alphaEminus)) - 1, min(min(alphaEminus))])

indexVecHelp = (1 : K).';

for nn = 1 : 3
  for np = 1 : 3
    indexVec = markE0TE0T{nn, np} * indexVecHelp;
    helperVec = (indexVec > 0);
    alphaEface(helperVec,nn) = min( alphaEplus(helperVec), alphaEminus(indexVec(helperVec)) ) .* (intFe(helperVec,nn) < 0) ...
        + min( alphaEminus(helperVec), alphaEplus(indexVec(helperVec)) ) .* (intFe(helperVec,nn) > 0) ;
  end % for np
end % for nn

% smoother_fun = @(x) 2*x.^4 - 5*x.^3 + 3*x.^2 + x;
% smoother_fun = @(x) 0.5 * (1-cos(x/pi));
% smoother_fun = @(x) -2 * x.^3 + 3 * x.^2;
% smoother_fun = @(x) x.^3;
smoother_fun = @(x) x;
alphaEface = smoother_fun(alphaEface);

assert(all(all((0 <= alphaEface) & (alphaEface <= 1))), 'alphaEface: %.3e %.3e', [max(max(alphaEface)) - 1, min(min(alphaEface))])

end % function getFluxLimitingParameter