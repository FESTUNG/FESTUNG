% Second step of the four-part algorithm in the main loop.

%===============================================================================
%> @file template/solveStep.m
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

% Obtain Runge-Kutta rule
[problemData.timeLvls, problemData.omega] = rungeKuttaExplicit(problemData.sweData.schemeOrder, problemData.sweData.dt, problemData.sweData.t);

% Initialize time stepping
K = problemData.sweData.K;
N = problemData.sweData.N;

problemData.sweData.tLvls = problemData.timeLvls;
problemData.sweData.omega = problemData.omega;

problemData.sweData.cDiscRK0 = [ reshape(problemData.sweData.cDisc(:,:,1).', K*N, 1) ; ...
                                 reshape(problemData.sweData.cDisc(:,:,2).', K*N, 1) ; ...
                                 reshape(problemData.sweData.cDisc(:,:,3).', K*N, 1) ];
problemData.sweData.cDiscRK = problemData.sweData.cDiscRK0;

problemData.transportData.timeLvls = problemData.timeLvls;
problemData.transportData.omega = problemData.omega;

problemData.transportData.cDiscRK0 = cellfun(@(c) reshape(c.', [K*N,1]), problemData.transportData.cDisc, 'UniformOutput', false);
problemData.transportData.cDiscRK = problemData.transportData.cDiscRK0;

% Perform RK stepping
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);
end % function