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
function pd = solveStep(pd, ~)
K = pd.K;
N = pd.N;
dt = pd.dt;
       
% Build system matrix and right hand side vectors
sysA = [ sparse(K*N,3*K*N); ...
         pd.tidalTerms{1}, sparse(K*N,2*K*N); ...
         pd.tidalTerms{2}, sparse(K*N,2*K*N) ];
sysV = [ pd.globL{1} - pd.globLRI{1}; ...
         pd.globL{2} - pd.globLRI{2}; ...
         pd.globL{3} - pd.globLRI{3} ] - ...
       [ sparse(K*N,1); ...
         pd.nonlinearTerms + pd.bottomFrictionTerms] - ... 
       pd.riemannTerms;

% Linearize solution vector
sysY = [ reshape(pd.cDisc(:,:,1).', K*N, 1) ; ...
         reshape(pd.cDisc(:,:,2).', K*N, 1) ; ...
         reshape(pd.cDisc(:,:,3).', K*N, 1) ];
       
% Compute solution at next time step using explicit or semi-implicit scheme
switch pd.scheme
  case 'explicit'
    sysY = sysY + dt * ( pd.sysW \ (sysV - (pd.linearTerms - sysA) * sysY) );

  case 'semi-implicit'
    sysY = (pd.sysW + dt * (pd.linearTerms - sysA)) \ (pd.sysW * sysY + dt * sysV);
          
  otherwise
    error('Invalid time-stepping scheme')
end % switch

% Reshape linearized vector to solution vectors
pd.cDisc(:,:,1) = reshape(sysY(        1 :   K*N), N, K).';
pd.cDisc(:,:,2) = reshape(sysY(  K*N + 1 : 2*K*N), N, K).';
pd.cDisc(:,:,3) = reshape(sysY(2*K*N + 1 : 3*K*N), N, K).';
end % function