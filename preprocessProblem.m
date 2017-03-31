% Performs all pre-processing steps, such as grid generation, assembly of
% stationary blocks, etc. for the problem solution.

%===============================================================================
%> @file darcyVert_sweVert/preprocessProblem.m
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
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
problemData.darcyData = problemData.darcySteps.preprocessProblem(problemData.darcyData);
problemData.sweData = problemData.sweSteps.preprocessProblem(problemData.sweData);

% Coupling matrices for hydraulic head
g1D = problemData.generateGrid1D(problemData.darcyData.numElem(1), problemData.darcyData.g);
problemData.darcyData.tildeGlobQcouple = assembleMatEdgeTetraPhiIntPhi1DIntNu(problemData.darcyData.g, g1D, problemData.darcyData.g.markE0TbdrCoupling, problemData.sweData.tildeHatQdiag);
problemData.darcyData.tildeGlobScouple = assembleMatEdgeTetraPhiIntPhi1DInt(problemData.darcyData.g, g1D, problemData.darcyData.g.markE0TbdrCoupling, problemData.sweData.tildeHatQdiag);

problemData.hatS = integrateRefEdgeTetraPhiIntPhiIntPhiExtPhiExt(problemData.sweData.N, problemData.sweData.qOrd, problemData.sweData.basesOnQuad2D);

% Grid data structure for coupling terms
problemData.gCoupling.areaE0T = problemData.sweData.g.areaE0T;
problemData.gCoupling.nuE0T = problemData.sweData.g.nuE0T;
% Edge 1 of SWE is coupled to edge 2 of Darcy - identified via the 1D-elements
problemData.gCoupling.markE0TE0T = cellfun(@(c) sparse(problemData.sweData.g.numT, problemData.darcyData.g.numT), cell(1,4), 'UniformOutput', false);
problemData.gCoupling.markE0TE0T{1} = bsxfun(@and, problemData.sweData.g.g1D.markT2DT, problemData.sweData.g.markE0TbdrCoupling(:, 1)) * ...
                                        bsxfun(@times, g1D.markT2DT, problemData.darcyData.g.markE0TbdrCoupling(:,2))';
problemData.gCoupling.markE0TE0T = cellfun(@(c) logical(c), problemData.gCoupling.markE0TE0T, 'UniformOutput', false);                                      
end