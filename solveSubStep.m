% Second step of the three-part algorithm in the iterateSubSteps loop of each
% step in the main loop.

%===============================================================================
%> @file template/solveSubStep.m
%>
%> @brief Second step of the three-part algorithm in the iterateSubSteps loop of
%>        each step in the main loop.
%===============================================================================
%>
%> @brief Second step of the three-part algorithm in the iterateSubSteps loop of
%>        each step in the main loop.
%>
%> The iterateSubSteps loop repeatedly executes three steps until the number of
%> substep iterations equals the order of the underlying Runge-Kutta method.
%> These three steps are:
%>
%>  1. preprocessSubStep()
%>  2. solveSubStep()
%>  3. postprocessSubStep()
%> 
%> This routine is executed second in each substep loop iteration and is
%> intended to produce the solution at the next substep, e.g., at the time-level
%> of the new Runge-Kutta step.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nSubStep     The current iteration number of the iterateSubSteps 
%>                      loop. 
%>
%> @retval problemData  The input struct enriched with solution data at
%>                      the next substep. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>                      Modified 08/17/16 by Hennes Hajduk
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
%
function problemData = solveSubStep(problemData, ~, nSubStep)
K = problemData.K;
N = problemData.N;
dt = problemData.dt;

% Compute height from potentially different approximation orders:
% zbDiscLin is always linear (N=3) while cDisc(:,:,1) can be of any
% approximation order
hDisc = reshape(computeOpDataDiscDataDisc(@minus, problemData.cDisc(:,:,1), problemData.zbDiscLin).', [], 1);

% Build right hand side vector
sysV = cell2mat(problemData.globL) - cell2mat(problemData.globLRI) - ...
       [ sparse(K*N,1); problemData.nonlinearTerms + problemData.bottomFrictionTerms] - ... 
       problemData.riemannTerms;

% Compute solution at next time step using explicit or semi-implicit scheme
switch problemData.schemeType
  case 'explicit'
    sysA = [ sparse(K*N,K*max(N,3)); problemData.tidalTerms{1}; problemData.tidalTerms{2} ];
    cDiscDot = problemData.sysW \ (sysV - problemData.linearTerms * problemData.cDiscRK + sysA * hDisc);
    
    % Apply slope limiting to time derivative
    for i = 1 : length(problemData.slopeLimList)
      switch problemData.slopeLimList{i}
        case 'elevation'
          cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot(1:K*N), [N K])', problemData.globM, problemData.globMDiscTaylor);
          cDiscDotTaylorLim = applySlopeLimiterTaylor(problemData.g, cDiscDotTaylor, problemData.g.markV0TbdrD, NaN(K,3), problemData.basesOnQuad, problemData.typeSlopeLim);
          cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + problemData.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
          cDiscDot(1:K*N) = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', problemData.globM, problemData.globMDiscTaylor)', [K*N 1]);
        case 'momentum'
          cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot(K*N+1:2*K*N), [N K])', problemData.globM, problemData.globMDiscTaylor);
          cDiscDotTaylorLim = applySlopeLimiterTaylor(problemData.g, cDiscDotTaylor, problemData.g.markV0TbdrD, NaN(K,3), problemData.basesOnQuad, problemData.typeSlopeLim);
          cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + problemData.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
          cDiscDot(K*N+1:2*K*N) = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', problemData.globM, problemData.globMDiscTaylor)', [K*N 1]);
          
          cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot(2*K*N+1:3*K*N), [N K])', problemData.globM, problemData.globMDiscTaylor);
          cDiscDotTaylorLim = applySlopeLimiterTaylor(problemData.g, cDiscDotTaylor, problemData.g.markV0TbdrD, NaN(K,3), problemData.basesOnQuad, problemData.typeSlopeLim);
          cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + problemData.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
          cDiscDot(2*K*N+1:3*K*N) = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', problemData.globM, problemData.globMDiscTaylor)', [K*N 1]);
        otherwise
          error('Slope limiting not implemented for non primary variables.')
      end % switch
    end % for
    
    problemData.cDiscRK = problemData.omega(nSubStep) * problemData.cDiscRK0 + (1 - problemData.omega(nSubStep)) * (problemData.cDiscRK + dt * cDiscDot);

  case 'semi-implicit'
    sysA = [ sparse(K*N,3*K*N); problemData.tidalTerms{1}, sparse(K*N,2*K*N); problemData.tidalTerms{2}, sparse(K*N,2*K*N) ];
    problemData.cDiscRK = (problemData.sysW + dt * (problemData.linearTerms - sysA)) \ (problemData.sysW * problemData.cDiscRK + dt * sysV);

  otherwise
    error('Invalid time-stepping scheme')
end % switch
end % function
