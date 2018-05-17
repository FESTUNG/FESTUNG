% Second step of the four-part algorithm in the main loop.

%===============================================================================
%> @file
%>
%> @brief Second step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until a certain
%> criterion, such as reaching the end of the simulation time or
%> steady state. These four steps are:
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
function problemData = solveStep(problemData, nStep)
K = problemData.K;
N = problemData.N;

if problemData.isAdaptiveTimestep
  problemData = selectTimeStep(problemData, nStep);
end % if

% Obtain Runge-Kutta rule
switch problemData.schemeType
  case 'explicit'
    [problemData.tLvls, problemData.omega] = rungeKuttaExplicit(problemData.schemeOrder, problemData.dt, problemData.t);

  case 'semi-implicit'
    problemData.tLvls = problemData.t + problemData.dt;
    
  otherwise
    error('Invalid time stepping scheme')
end % switch

% Initialize solution vectors for RK steps
problemData.cDiscRK0 = [reshape(problemData.cDisc(:,:,1).', K*N, 1); ...
                        reshape(problemData.cDisc(:,:,2).', K*N, 1); ...
                        reshape(problemData.cDisc(:,:,3).', K*N, 1)];
problemData.cDiscRK = problemData.cDiscRK0;

% Carry out RK steps
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);

end % function