function [highOrderSol, success] = Zalesak(highOrderSol, lowOrderMeans, suppressedFluxes, markE0TE0T, K, N, tau, areaT, maxCounter)

%% Set parameters.
Cmax = +1; Cmin = -1; epsilon = 1e-8; residualCrit = 1e-7;% maxCounter = 100;

%% Initialize alphaEface and First residual
alphaEface = zeros(K,3);
residual = norm(norm(suppressedFluxes, Inf), Inf);
counter = 0;

%% Start iteration
while residual > residualCrit && counter < maxCounter
  
  counter = counter + 1;
  oldSuppressedFluxes = suppressedFluxes;
  
  % Compute the alpha_E
  assert(all((Cmin - epsilon <= lowOrderMeans) & (lowOrderMeans <= Cmax + epsilon)), 'cMean: %.3e %.3e', [max(max(lowOrderMeans)), min(min(lowOrderMeans))])
  alphaEplus = min( 1, areaT .* max(0, Cmax - lowOrderMeans) ./ ( -sum( suppressedFluxes .* (suppressedFluxes < 0) , 2 ) * tau + epsilon) );
  alphaEminus = min( 1, areaT .* max(0, lowOrderMeans - Cmin) ./ ( sum( suppressedFluxes .* (suppressedFluxes > 0) , 2 ) * tau + epsilon) );
  
  % Compute the alpah_E,E'
  indexVecHelp = (1 : K).';
  for nn = 1 : 3
    for np = 1 : 3
      indexVec = markE0TE0T{nn, np} * indexVecHelp;
      helperVec = (indexVec > 0);
      alphaEface(helperVec,nn) = min( alphaEplus(helperVec), alphaEminus(indexVec(helperVec)) ) .* (suppressedFluxes(helperVec,nn) < 0) ...
        + min( alphaEminus(helperVec), alphaEplus(indexVec(helperVec)) ) .* (suppressedFluxes(helperVec,nn) > 0) ;
    end % for np
  end % for nn
  
  assert(all(all((0 <= alphaEface) & (alphaEface <= 1))), 'alphaEface: %.3e %.3e', [max(max(alphaEface)) - Cmax, min(min(alphaEface)) - Cmin])
  
  % Compute the updated solution and fluxes
  for n = 1 : 3
    lowOrderMeans = lowOrderMeans - tau * alphaEface(:,n) .* suppressedFluxes(:,n) ./ areaT;
    suppressedFluxes(:,n) = (1 - alphaEface(:,n)) .* suppressedFluxes(:,n);
  end
  
  % Check new residual
  residual = norm(norm(suppressedFluxes, Inf), Inf);
  deltaSuppressed = norm(norm(oldSuppressedFluxes - suppressedFluxes, Inf), Inf);
  if deltaSuppressed < residualCrit
    break;
  end
  
end % while norm > residualCrit

highOrderSol(1:N:K*N) = lowOrderMeans ./ sqrt(2);

success = (residual < residualCrit);

fprintf('Zalesak: residual = %d, counter = %d, delta suppressed flux = %d.\n', residual, counter, deltaSuppressed);

end % function getFluxLimitingParameter