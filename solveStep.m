% Second step of the four-part algorithm in the main loop. Computes the
% solution at the next time step.

%===============================================================================
%> @file diffusion/solveStep.m
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>        Computes the solution at the next time step.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>        Computes the solution at the next time step.
%>
%> The main loop repeatedly executes four steps until the parameter
%> <code>problemData.isFinished</code> becomes <code>true</code>.
%> These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed second in each loop iteration.
%> It assembles the global system and solves it using the backslash-operator.
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
function problemData = solveStep(problemData, nStep) %#ok<INUSD>
K = problemData.K;
N = problemData.N;
%% Building and solving the system.
sysA = [                                                problemData.globM,                                                  sparse(K*N,K*N), -problemData.globH{1}+problemData.globQ{1}+problemData.globQN{1};
                                                          sparse(K*N,K*N),                                                problemData.globM, -problemData.globH{2}+problemData.globQ{2}+problemData.globQN{2};
         -problemData.globG{1}+problemData.globR{1}+problemData.globRD{1}, -problemData.globG{2}+problemData.globR{2}+problemData.globRD{2},                             problemData.globS+problemData.globSD];
sysV = [-problemData.globJD{1}; -problemData.globJD{2}; problemData.globKD-problemData.globKN+problemData.globL];
problemData.sysY = (problemData.sysW + problemData.tau*sysA) \ (problemData.sysW*problemData.sysY + problemData.tau*sysV);
end