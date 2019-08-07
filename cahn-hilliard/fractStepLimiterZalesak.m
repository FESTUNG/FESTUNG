% Fractional step flux limiter to correct local mean values using Zalesak's
% flux-corrected transport algorithm for the computation of worst-case
% correction factors.

%===============================================================================
%> @file
%>
%> @brief % Fractional step flux limiter to correct local mean values using Zalesak's
%>          flux-corrected transport algorithm for the computation of worst-case
%>          correction factors.
%===============================================================================
%>
%> @brief % Fractional step flux limiter to correct local mean values using Zalesak's
%>          flux-corrected transport algorithm for the computation of worst-case
%>          correction factors.
%>
%> Foo bar
%>
%> @param  g             Grid/triangulation (see 
%>                       <code>generateGridData()</code>) 
%>                       @f$[1 \times 1 \text{ struct}]@f$
%> @param  umin          Lower solution bound.
%> @param  umax          Upper solution bound.
%> @param  highOrderSol  DG/modal basis representation of high-order solution. @f$[K\cdot N \times 1]@f$
%> @param  lowOrderMeans Mean values of low-order solution. @f$[K \times 1]@f$
%>
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2019 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Andreas Rupp, Florian Frank, 2019
%> 
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock

function [highOrderSol, success] = fractStepLimiterZalesak(g, umin, umax, highOrderSol, lowOrderMeans, ...
  suppressedFluxes, K, N, tau, maxIter)

%% Asserts.
assert(umin < umax)
assert(isequal(size(highOrderSol), [K*N, 1]))
assert(isequal(size(lowOrderMeans), [K, 1]))

%% Set parameters.
epsilon = 1e-8; 
epsilonSupprDelta = 1e-7;
epsilonSupprFlux  = 1e-7;

%% Initialize alphaEface and initial residual.
alphaEface = zeros(K, 3);
normSupprFlux = norm(norm(suppressedFluxes, Inf), Inf);
numIter = 0;

%% Start iteration.
while (normSupprFlux > epsilonSupprFlux) && (numIter < maxIter)
  
  numIter = numIter + 1;
  oldSuppressedFluxes = suppressedFluxes;
  
  % Compute element correction factors alpha_E.
  assert(all((umin - epsilon <= lowOrderMeans) & (lowOrderMeans <= umax + epsilon)), 'uMean: %.3e %.3e', [max(max(lowOrderMeans)), min(min(lowOrderMeans))])
  alphaEplus = min(1, g.areaT .* max(0, umax - lowOrderMeans) ./ ( -sum( suppressedFluxes .* (suppressedFluxes < 0) , 2 ) * tau + epsilon) );
  alphaEminus = min(1, g.areaT .* max(0, lowOrderMeans - umin) ./ ( sum( suppressedFluxes .* (suppressedFluxes > 0) , 2 ) * tau + epsilon) );
  
  % Compute edge correction factors alpha_E,E'.
  indexVecHelp = (1 : K).';
  for nn = 1 : 3
    for np = 1 : 3
      indexVec = g.markE0TE0T{nn, np} * indexVecHelp;
      helperVec = (indexVec > 0);
      alphaEface(helperVec,nn) = min( alphaEplus(helperVec), alphaEminus(indexVec(helperVec)) ) .* (suppressedFluxes(helperVec,nn) < 0) ...
        + min( alphaEminus(helperVec), alphaEplus(indexVec(helperVec)) ) .* (suppressedFluxes(helperVec,nn) > 0) ;
    end % for np
  end % for nn
  
  assert(all(all((0 <= alphaEface) & (alphaEface <= 1))), 'alphaEface: %.3e %.3e', [max(max(alphaEface)) - umax, min(min(alphaEface)) - umin])
  
  % Compute the updated solution and fluxes.
  for n = 1 : 3
    lowOrderMeans = lowOrderMeans - tau * alphaEface(:,n) .* suppressedFluxes(:,n) ./ g.areaT;
    suppressedFluxes(:,n) = (1 - alphaEface(:,n)) .* suppressedFluxes(:,n);
  end
  
  % Check new normSupprFlux.
  normSupprFlux  = norm(norm(suppressedFluxes, Inf), Inf);
  normSupprDelta = norm(norm(oldSuppressedFluxes - suppressedFluxes, Inf), Inf);
  if normSupprDelta < epsilonSupprDelta
    break
  end
  
end % while norm > epsilonSupprFlux

highOrderSol(1:N:K*N) = lowOrderMeans ./ sqrt(2); % sqrt(2) is the value of the first trial function

success = (normSupprFlux < epsilonSupprFlux);

% Todo: The delta should only be printed if suppressed flux is unequal to
% zero.
fprintf('Zalesak: suppressed flux: %d, iterations: %d, delta suppressed flux: %d.\n', normSupprFlux, numIter, normSupprDelta);

end % function fractStepLimiterZalesak
