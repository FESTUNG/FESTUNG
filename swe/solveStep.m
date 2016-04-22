% Second step of the four-part algorithm in the main loop.

%===============================================================================
%> @file template/solveStep.m
%>
%> @brief Second step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until the number of
%> iterations provided by configureProblem in the parameter
%> <code>numSteps</code> is reached. These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed second in each loop iteration and is intended to
%> produce the solution at the next step, e.g., at a new time-level.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with solution data at
%>                      the next step. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function problemData = solveStep(problemData, ~)
K = problemData.K;
N = problemData.N;
dt = problemData.dt;
       
% Build system matrix and right hand side vectors
if problemData.isOSRiem
  sysV = [ problemData.globL{1} + problemData.globVOSRiem{2}; ...
           problemData.globL{2} + 0.5 * problemData.globROSold{1}; ...
           problemData.globL{3} + 0.5 * problemData.globROSold{2} ];
  sysA = [ problemData.globVOSRiem{1}, sparse(K*N,2*K*N) ; ...
                      sparse(K*N,K*N), 0.5 * problemData.globUOSold,              sparse(K*N,K*N) ; ...
                      sparse(K*N,K*N),              sparse(K*N,K*N), 0.5 * problemData.globUOSold ];    
else
  sysV = [ problemData.globL{1}; ...
           problemData.globL{2} + problemData.globROSold{1}; ...
           problemData.globL{3} + problemData.globROSold{2} ];
  sysA = [ sparse(K*N,3*K*N) ; ...
           sparse(K*N,K*N), problemData.globUOSold, sparse(K*N,K*N); ...
           sparse(K*N,K*N),        sparse(K*N,K*N),   globUOSold ];
end % if

% Linearize solution vector
sysY = [ reshape(problemData.cDisc(:,:,1).', K*N, 1) ; ...
         reshape(problemData.cDisc(:,:,2).', K*N, 1) ; ...
         reshape(problemData.cDisc(:,:,3).', K*N, 1) ];
       
% Compute solution at next time step using explicit or semi-implicit scheme
switch problemData.scheme
  case 'explicit'
    sysY = sysY + dt * ( problemData.sysW \ (sysV - ((sysA + problemData.linearTerms) * sysY + ...
            [sparse(K*N,1); problemData.nonlinearTerms + problemData.bottomFrictionTerms] + ...
            problemData.riemannTerms) ) );
  case 'semi-implicit'
    sysY = (problemData.sysW + dt * problemData.linearTerms) \ ( dt * sysV - dt * (sysA * sysY + ...
            [sparse(K*N,1); problemData.nonlinearTerms + problemData.bottomFrictionTerms] + ...
            problemData.riemannTerms) + problemData.sysW * sysY );
  otherwise
    error('Invalid time-stepping scheme')
end % switch

% Reshape linearized vector to solution vectors
problemData.cDisc(:,:,1) = reshape(sysY(        1 :   K*N), N, K).';
problemData.cDisc(:,:,2) = reshape(sysY(  K*N + 1 : 2*K*N), N, K).';
problemData.cDisc(:,:,3) = reshape(sysY(2*K*N + 1 : 3*K*N), N, K).';
end % function