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
function pd = solveSubStep(pd, ~, nSubStep)
K = pd.K;
N = pd.N;
dt = pd.dt;

% Compute height from potentially different approximation orders:
% zbDiscLin is always linear (N=3) while cDisc(:,:,1) can be of any
% approximation order
hDisc = reshape(execin('swe/computeSumDataDiscDataDisc', pd.cDisc(:,:,1), -pd.zbDiscLin).', [], 1);

% Build right hand side vector
sysV = cell2mat(pd.globL) - cell2mat(pd.globLRI) - ...
       [ sparse(K*N,1); pd.nonlinearTerms + pd.bottomFrictionTerms] - ... 
       pd.riemannTerms;

% Compute solution at next time step using explicit or semi-implicit scheme
switch pd.schemeType
  case 'explicit'
    sysA = [ sparse(K*N,K*max(N,3)); pd.tidalTerms{1}; pd.tidalTerms{2} ];
    cDiscDot = pd.sysW \ (sysV - pd.linearTerms * pd.cDiscRK + sysA * hDisc );
    pd.cDiscRK = pd.omega(nSubStep) * pd.cDiscRK0 + (1 - pd.omega(nSubStep)) * (pd.cDiscRK + dt * cDiscDot);

  case 'semi-implicit'
    sysA = [ sparse(K*N,3*K*N); pd.tidalTerms{1}, sparse(K*N,2*K*N); pd.tidalTerms{2}, sparse(K*N,2*K*N) ];
    pd.cDiscRK = (pd.sysW + dt * (pd.linearTerms - sysA)) \ (pd.sysW * pd.cDiscRK + dt * sysV);
          
  otherwise
    error('Invalid time-stepping scheme')
end % switch
end % function
