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
%> @retval highOrderSol  Flux-corrected high-order solution (in modal basis) @f$[K\cdot N \times 1]@f$.
%> @retval isSuccess
%> @param  g             Grid/triangulation (see 
%>                       <code>generateGridData()</code>) 
%>                       @f$[1 \times 1 \text{ struct}]@f$
%> @param  umin          Lower solution bound.
%> @param  umax          Upper solution bound.
%> @param  highOrderSol  DG/modal basis representation of high-order solution @f$[K\cdot N \times 1]@f$.
%> @param  lowOrderMeans Mean values of low-order solution. @f$[K \times 1]@f$
%> @param  supprFluxes   Target flux, i.e., difference of high order and low order fluxes [???].
%> @param  N             Local degrees of freedom.
%> @param  epsSupprFlux  Tolerance for the suppressed flux.
%> @param  epsSupprDelta Tolerance for the abs. distance of suppressed fluxes in two iterations.
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

function [highOrderSol, isSuccess] = fractStepLimiterZalesak(g, umin, umax, ...
  highOrderSol, lowOrderMeans, supprFluxes, N, tau, ...
  epsSupprFlux, epsSupprDelta, maxIter)

%% Set parameters.
deltaBound = 1e-8; 
deltaDenom = 1e-8;
K = g.numT;

%% Asserts.
assert(umin < umax)
assert(isequal(size(highOrderSol), [K*N, 1]))
assert(isequal(size(lowOrderMeans), [K, 1]))

%% Initialize alphaEface and initial residual.
alphaEface = zeros(K, 3);
normSupprFlux = norm(norm(supprFluxes, Inf), Inf);
numIter = 0;

%% Start iteration.
while (normSupprFlux > epsSupprFlux) && (numIter < maxIter)
  
  numIter = numIter + 1;
  oldSuppressedFluxes = supprFluxes;
  
  % Verify that all means are within [umin, umax] (up to arithmetic errors).
  assert(all((umin - deltaBound <= lowOrderMeans) & (lowOrderMeans <= umax + deltaBound)), ...
    'Mean values exceed bounds (possibly in small digits): %.3e %.3e', max(max(lowOrderMeans)), min(min(lowOrderMeans)))
  
  % Compute element correction factors alpha_E.
  alphaEplus  = min(1, g.areaT .* max(0, umax - lowOrderMeans) ./ (-sum( supprFluxes .* (supprFluxes < 0), 2 ) * tau + deltaDenom) );
  alphaEminus = min(1, g.areaT .* max(0, lowOrderMeans - umin) ./ ( sum( supprFluxes .* (supprFluxes > 0), 2 ) * tau + deltaDenom) );
  
  % Compute edge correction factors alpha_E,E'.
  indexVecHelp = (1 : K).';
  for nn = 1 : 3
    for np = 1 : 3
      indexVec = g.markE0TE0T{nn, np} * indexVecHelp;
      helperVec = (indexVec > 0);
      alphaEface(helperVec,nn) = min( alphaEplus(helperVec), alphaEminus(indexVec(helperVec)) ) .* (supprFluxes(helperVec,nn) < 0) ...
        + min( alphaEminus(helperVec), alphaEplus(indexVec(helperVec)) ) .* (supprFluxes(helperVec,nn) > 0) ;
    end % for np
  end % for nn
  
  % Verify that all correction factors are within [0,1].
  assert(all(all((0 <= alphaEface) & (alphaEface <= 1))), 'alphaEface: %.3e %.3e', max(max(alphaEface)) - umax, min(min(alphaEface)) - umin)
  
  % Compute the updated solution and fluxes.
  for n = 1 : 3
    lowOrderMeans = lowOrderMeans - tau * alphaEface(:,n) .* supprFluxes(:,n) ./ g.areaT;
    supprFluxes(:,n) = (1 - alphaEface(:,n)) .* supprFluxes(:,n);
  end
  
  % Check new normSupprFlux.
  normSupprFlux  = norm(norm(supprFluxes, Inf), Inf);
  normSupprDelta = norm(norm(oldSuppressedFluxes - supprFluxes, Inf), Inf);
  if normSupprDelta < epsSupprDelta
    break
  end
  
end % while norm > epsSupprFlux

highOrderSol(1:N:K*N) = lowOrderMeans ./ sqrt(2); % sqrt(2) is the value of the first trial function
isSuccess = (normSupprFlux < epsSupprFlux);

% Todo: The delta should only be printed if suppressed flux is unequal to
% zero.
fprintf('Zalesak: suppressed flux: %d, iterations: %d, delta suppressed flux: %d.\n', normSupprFlux, numIter, normSupprDelta);

end % function fractStepLimiterZalesak
