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
%> The matrices @f$\mathsf{{S}}^n \in \mathbb{R}^{KN\times KR}@f$are defined as 
%> @f[
%> [\mathsf{{S}}^n]_{(k-1)N+i,(k-1)N+r} =
%>  w_r \varphi_{ki}(\mathbf{q}_r) \,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times n_\mathrm{edges}]@f$
%> @param refEdgePhiIntPerQuad  Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPerQuad()</code>.
%>                    @f$[N \times R \times n_\mathrm{edges}]@f$
%> @retval ret        The assembled matrices @f$[ n_\mathrm{edges}\times1 \mathrm{cell}]@f$
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

function [highOrderSol, success] = fractStepLimiterZalesak(highOrderSol, lowOrderMeans, suppressedFluxes, markE0TE0T, K, N, tau, areaT, maxCounter)

% K, N, tau, areaT are elements of problemData

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

highOrderSol(1:N:K*N) = lowOrderMeans ./ sqrt(2); % sqrt(2) is the value of the first trial function

success = (residual < residualCrit);

fprintf('Zalesak: residual = %d, counter = %d, delta suppressed flux = %d.\n', residual, counter, deltaSuppressed);

end % function getFluxLimitingParameter
