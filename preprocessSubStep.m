% First step of the three-part algorithm in the iterateSubSteps loop of each
% step in the main loop.
%
%===============================================================================
%> @file template/preprocessSubStep.m
%>
%> @brief First step of the three-part algorithm in the iterateSubSteps loop of 
%>        each step in the main loop.
%===============================================================================
%>
%> @brief First step of the three-part algorithm in the iterateSubSteps loop of 
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
%> This routine is executed first in each substep loop iteration and is intended
%> to execute preprocessing operations, e.g., evaluate boundary conditions or
%> right hand side values, assemble time-dependent matrix blocks, etc.
%>
%> @param  pd           A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the iterateSubSteps 
%>                      loop. 
%>
%> @retval pd           The input struct enriched with preprocessed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
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
function problemData = preprocessSubStep(problemData, nStep, nSubStep)
problemData.sweData = problemData.swe_preprocessSubStep(problemData.sweData, nStep, nSubStep);

% use velocities of swe for transport problem
problemData.transportData.vNormalOnQuadEdge = problemData.sweData.massFluxQ0E0T;
problemData.transportData.hDisc = problemData.sweData.hDisc;
problemData.transportData.uHDisc = problemData.sweData.uHDisc;
problemData.transportData.vHDisc = problemData.sweData.vHDisc;

problemData.transportData = problemData.transport_preprocessSubStep(problemData.transportData, nStep, nSubStep);
end % function
