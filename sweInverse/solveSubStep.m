% Second step of the three-part algorithm in the iterateSubSteps loop of each
% step in the main loop.

%===============================================================================
%> @file sweInverse/solveSubStep.m
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
%> @param  pd           A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nSubStep     The current iteration number of the iterateSubSteps 
%>                      loop. 
%>
%> @retval pd           The input struct enriched with solution data at
%>                      the next substep. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function pd = solveSubStep(pd, ~, nSubStep)
K = pd.K;
N = pd.N;
dt = pd.dt;

diffusionTerms = sparse(3*K*N,1);
if pd.isDiffusion
  for m = 1:2
    q = pd.sysW \ ( pd.globDiffusion{m} * pd.cDiscRK - [pd.rhsDiffusion{1}{m}; pd.rhsDiffusion{2}{m}; pd.rhsDiffusion{3}{m}] );
    diffusionTerms = diffusionTerms - pd.globDiffusion2Mom{m} * q;
  end % for
end % if

% Build right hand side vector
sysV = cell2mat(pd.globL) - cell2mat(pd.globLRI) - cell2mat(pd.globLF) - ...
       [ -pd.nudgingTerm + pd.penalty * reshape(pd.zbDisc.', N*K, 1); pd.nonlinearTerms + pd.bottomFrictionTerms ] ... 
       - pd.riemannTerms - diffusionTerms;

sysA = pd.gravityTerm - [pd.tidalTerms{1}; pd.tidalTerms{2}];

% Compute solution at next time step using explicit or semi-implicit scheme
switch pd.schemeType
  case 'explicit'
    cDiscDot = pd.sysW \ (sysV - [ [sparse(K*N,K*N); sysA], pd.linearTerms ] * pd.cDiscRK);
    
    % Apply slope limiting to time derivative
    for i = 1 : length(pd.slopeLimList)
      switch pd.slopeLimList{i}
        case 'elevation'
          cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot(1:K*N), [N K])', pd.globM, pd.globMDiscTaylor);
          cDiscDotTaylorLim = applySlopeLimiterTaylor(pd.g, cDiscDotTaylor, pd.g.markV0TbdrD, NaN(K,3), pd.basesOnQuad, pd.typeSlopeLim);
          cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + pd.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
          cDiscDot(1:K*N) = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', pd.globM, pd.globMDiscTaylor)', [K*N 1]);
        case 'momentum'
          cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot(K*N+1:2*K*N), [N K])', pd.globM, pd.globMDiscTaylor);
          cDiscDotTaylorLim = applySlopeLimiterTaylor(pd.g, cDiscDotTaylor, pd.g.markV0TbdrD, NaN(K,3), pd.basesOnQuad, pd.typeSlopeLim);
          cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + pd.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
          cDiscDot(K*N+1:2*K*N) = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', pd.globM, pd.globMDiscTaylor)', [K*N 1]);
          
          cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot(2*K*N+1:3*K*N), [N K])', pd.globM, pd.globMDiscTaylor);
          cDiscDotTaylorLim = applySlopeLimiterTaylor(pd.g, cDiscDotTaylor, pd.g.markV0TbdrD, NaN(K,3), pd.basesOnQuad, pd.typeSlopeLim);
          cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + pd.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
          cDiscDot(2*K*N+1:3*K*N) = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', pd.globM, pd.globMDiscTaylor)', [K*N 1]);
        otherwise
          error('Slope limiting not implemented for non primary variables.')
      end % switch
    end % for
    
    pd.cDiscRK = pd.omega(nSubStep) * pd.cDiscRK0 + (1 - pd.omega(nSubStep)) * (pd.cDiscRK + dt * cDiscDot);

  case 'semi-implicit'
    sysA = [ sparse(K*N,3*K*N); pd.tidalTerms{1}, sparse(K*N,2*K*N); pd.tidalTerms{2}, sparse(K*N,2*K*N) ];
    pd.cDiscRK = (pd.sysW + dt * (pd.linearTerms + sysA)) \ (pd.sysW * pd.cDiscRK + dt * sysV);
  otherwise
    error('Invalid time-stepping scheme')
end % switch
end % function
