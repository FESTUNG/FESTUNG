% Fractional step flux limiter to correct local mean values using Zalesak's
% flux-corrected transport algorithm for the computation of worst-case
% correction factors.

%===============================================================================
%> @file
%>
%> @brief Fractional step flux limiter to correct local mean values using Zalesak's
%>        flux-corrected transport algorithm for the computation of worst-case
%>        correction factors.
%===============================================================================
%>
%> @brief Fractional step flux limiter to correct local mean values using Zalesak's
%>        flux-corrected transport algorithm for the computation of worst-case
%>        correction factors.
%>
%> The Algorithm is taken from Frank, Rupp, Kuzmin, 2019, "Bound-preserving flux 
%> limiting schemes for DG~discretizations of conservation laws with applications
%> to the Cahn-Hilliard equation".
%>
%> The implemented stopping criteria are:
%> - Criterion 0: @f$\max_{E,E'}|\mathcal{H}_{E,E'}^{(m)}| < \mathrm{epsSupprFlux}@f$
%> - Criterion 1: @f$\max_{E,E'}|\mathcal{H}_{E,E'}^{(m)} - \mathcal{H}_{E,E'}^{(m-1)}| < \mathrm{epsSupprDelta}@f$ 
%> - Criterion 2: @f$m = \mathrm{maxIter}@f$ 
%>
%>
%> @retval highOrderSol  Flux-corrected high-order solution (in modal basis) @f$[K\cdot N \times 1]@f$.
%> @retval flag          The number of the stopping criterion that triggered, see description.
%>
%> @param  g             Grid/triangulation (see 
%>                       <code>generateGridData()</code>) 
%>                       @f$[1 \times 1 \text{ struct}]@f$
%> @param  umin          Lower solution bound.
%> @param  umax          Upper solution bound.
%> @param  highOrderSol  DG/modal basis representation of high-order solution @f$[K\cdot N \times 1]@f$.
%> @param  lowOrderMeans Mean values of low-order solution. @f$[K \times 1]@f$
%> @param  supprFluxes   Target flux, i.e., difference of high order and low order fluxes @f$[K \times 3]@f$.
%> @param  N             Local degrees of freedom.
%> @param  tau           Time step size.
%> @param  epsSupprFlux  Tolerance for the suppressed flux.
%> @param  epsSupprDelta Tolerance for the abs. distance of suppressed fluxes in two iterations.
%> @param  maxIter       Maximum number of iterations; if reached,
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

function [highOrderSol, flag] = fractStepLimiterZalesak(g, umin, umax, ...
  highOrderSol, lowOrderMeans, supprFluxes, N, tau, ...
  epsSupprFlux, epsSupprDelta, maxIter)

%% Set parameters.
deltaBound = 10*eps(umin); % Arithmetic perturbations may violate the bounds, therefore, damping. 
K = g.numT;

%% Asserts.
assert(umin < umax)
assert(isequal(size(highOrderSol), [K*N, 1]))
assert(isequal(size(lowOrderMeans), [K, 1]))

%% Initialize alphaEface and initial residual.
alphaEface = zeros(K, 3);
normSupprFlux = norm(norm(supprFluxes, Inf), Inf);

%% Start iteration.
for numIter = 1 : maxIter
 
  % Save suppressed fluxes.
  oldSuppressedFluxes = supprFluxes;
  
  % Verify that all means are within [umin, umax].
  assert(all((umin <= lowOrderMeans) & (lowOrderMeans <= umax)), ...
    'Mean values exceed bounds before iteration %d (possibly in small digits): %.16f %.16f', numIter, min(min(lowOrderMeans)), max(max(lowOrderMeans)))
  
  % Compute element correction factors alpha_E. In Matlab, zero instead of eps 
  % would also work, since the denominator then is -0 + +0 = +0.
  alphaEplus  = min(1, g.areaT .* max(0, (umax - deltaBound) - lowOrderMeans) ./ (-sum( supprFluxes .* (supprFluxes < 0), 2) * tau + eps) );
  alphaEminus = min(1, g.areaT .* max(0, lowOrderMeans - (umin + deltaBound)) ./ ( sum( supprFluxes .* (supprFluxes > 0), 2) * tau + eps) );
  
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
  assert(all(all((0 <= alphaEface) & (alphaEface <= 1))), 'alphaEface: %.16f %.16f', min(min(alphaEface)), max(max(alphaEface)))
  
  % Compute the updated solution and fluxes.
  for n = 1 : 3
    lowOrderMeans = lowOrderMeans - tau * alphaEface(:,n) .* supprFluxes(:,n) ./ g.areaT;
    supprFluxes(:,n) = (1 - alphaEface(:,n)) .* supprFluxes(:,n);
  end
  
  % Compute new errors.
  normSupprFlux  = norm(norm(supprFluxes, Inf), Inf);
  normSupprDelta = norm(norm(oldSuppressedFluxes - supprFluxes, Inf), Inf);
  
  % Check stopping criteria.
  if normSupprFlux < epsSupprFlux
    flag = 0;
    break
  elseif normSupprDelta < epsSupprDelta
    flag = 1;
    break
  elseif numIter == maxIter
    flag = 2;
    break
  end
end % for

% Verify that all means are within [umin, umax].
assert(all((umin <= lowOrderMeans) & (lowOrderMeans <= umax)), ...
  'Mean values exceed bounds after iteration %d (possibly in small digits): %.16f %.16f', numIter, min(min(lowOrderMeans)), max(max(lowOrderMeans)))

% Scale coefficient vector.
highOrderSol(1:N:K*N) = lowOrderMeans ./ sqrt(2); % sqrt(2) is the value of the first trial function, i.e., phi(1,0,0).
% highOrderSol(1:N:K*N) = lowOrderMeans; 

% Print information.
fprintf('FSL: flag: %d, suppressed flux: %.3e, delta suppressed flux: %.3e, iterations: %d.\n', ...
  flag, normSupprFlux, normSupprDelta, numIter);

end % function fractStepLimiterZalesak
