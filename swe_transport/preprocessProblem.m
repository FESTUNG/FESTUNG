% Performs all pre-processing steps, such as grid generation, assembly of
% stationary blocks, etc. for the problem solution.

%===============================================================================
%> @file template/preprocessProblem.m
%>
%> @brief Performs all pre-processing tasks, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%===============================================================================
%>
%> @brief Performs all pre-processing steps, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%>
%> This routine is called after template/configureProblem.m.
%>
%> Typically, this step consists of grid generation, computation of derived
%> data structures, pre-computing often needed values (e.g., basis
%> functions on quadrature points), or assembly of time-independent matrix
%> blocks.
%>
%> @param  problemData  A struct with problem parameters, as provided by
%>                      configureProblem(). @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description and precomputed fields.
%>                      @f$[\text{struct}]@f$
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
function problemData = preprocessProblem(problemData)
h = getFunctionHandle('swe/preprocessProblem');
problemData.sweData = h(problemData.sweData);

problemData.transportData.K = problemData.sweData.K;
problemData.transportData.tau = problemData.sweData.dt;
problemData.transportData.velN = problemData.sweData.N;

h = getFunctionHandle('transport/preprocessProblem');
problemData.transportData = h(problemData.transportData);

% only created function handles for routines that are called repeatedly
problemData.swe_preprocessStep = getFunctionHandle('swe/preprocessStep');
problemData.swe_solveStep = getFunctionHandle('swe/solveStep');
problemData.swe_preprocessSubStep = getFunctionHandle('swe/preprocessSubStep');
problemData.swe_solveSubStep = getFunctionHandle('swe/solveSubStep');
problemData.swe_postprocessSubStep = getFunctionHandle('swe/postprocessSubStep');
problemData.swe_postprocessStep = getFunctionHandle('swe/postprocessStep');
problemData.swe_outputStep = getFunctionHandle('swe/outputStep');
problemData.transport_preprocessStep = getFunctionHandle('transport/preprocessStep');
problemData.transport_solveStep = getFunctionHandle('transport/solveStep');
problemData.transport_preprocessSubStep = getFunctionHandle('transport/preprocessSubStep');
problemData.transport_solveSubStep = getFunctionHandle('transport/solveSubStep');
problemData.transport_postprocessSubStep = getFunctionHandle('transport/postprocessSubStep');
problemData.transport_postprocessStep = getFunctionHandle('transport/postprocessStep');
problemData.transport_outputStep = getFunctionHandle('transport/outputStep');
end % function