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
function pd = solveStep(pd, nStep)
K = pd.K;
N = pd.N;

% Obtain Runge-Kutta rule
switch pd.schemeType
  case 'explicit'
    [pd.tLvls, pd.omega] = rungeKuttaSSP(pd.schemeOrder, pd.dt, pd.t);

  case 'semi-implicit'
    pd.tLvls = pd.t + dt;
    
  otherwise
    error('Invalid time stepping scheme')
end % switch
    
% Initialize solution vectors for RK steps
pd.cDiscRK0 = [ reshape(pd.cDisc(:,:,1).', K*N, 1) ; ...
                reshape(pd.cDisc(:,:,2).', K*N, 1) ; ...
                reshape(pd.cDisc(:,:,3).', K*N, 1) ];
pd.cDiscRK = pd.cDiscRK0;

% Compute height
pd.sysH = pd.cDiscRK0(1:K*N) - reshape(pd.zbDisc.', K*N,1);

% Carry out RK steps
pd.isSubSteppingFinished = false;
pd = iterateSubSteps(pd, nStep);

% Compute change
if pd.isSteadyState
  pd.changeL2 = norm(pd.cDiscRK - pd.cDiscRK0, 2);
end % if

% Reshape linearized vector to solution vectors
pd.cDisc(:,:,1) = reshape(pd.cDiscRK(        1 :   K*N), N, K).';
pd.cDisc(:,:,2) = reshape(pd.cDiscRK(  K*N + 1 : 2*K*N), N, K).';
pd.cDisc(:,:,3) = reshape(pd.cDiscRK(2*K*N + 1 : 3*K*N), N, K).';
end % function