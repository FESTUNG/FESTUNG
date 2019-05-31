% Calls the preprocessing routines of the subproblems and initializes grid
% data structures for the coupling.

%===============================================================================
%> @file
%>
%> @brief Calls the preprocessing routines of the subproblems and initializes 
%>        grid data structures for the coupling.
%===============================================================================
%>
%> @brief Calls the preprocessing routines of the subproblems and initializes 
%>        grid data structures for the coupling.
%>
%> This routine is called after darcy_swe_2dv/configureProblem.m .
%>
%> It executes darcy_2dv/preprocessProblem.m and @link
%> swe_2dv/preprocessProblem.m @endlink.
%> Furthermore, grid data structures and static matrices for the coupling
%> are initialized.
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
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018
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
problemData.darcyData.globVeeQcouple = assembleMatEdgeQuadriPhiIntPhi1DIntNu(problemData.darcyData.g, g1D, problemData.darcyData.g.markE0TbdrCoupling, problemData.sweData.hatVeeQdiag);
problemData.darcyData.globVeeScouple = assembleMatEdgeQuadriPhiIntPhi1DInt(problemData.darcyData.g, g1D, problemData.darcyData.g.markE0TbdrCoupling, problemData.sweData.hatVeeQdiag, ones(problemData.darcyData.g.numT, 4));

% Grid data structure for coupling terms
problemData.markE0TE0T = bsxfun(@and, problemData.sweData.g.g1D.markT2DT, problemData.sweData.g.markE0TbdrCoupling(:, 1)) * ...
                                        bsxfun(@times, g1D.markT2DT, problemData.darcyData.g.markE0TbdrCoupling(:,2))';
problemData.markT2DT = g1D.markT2DT;                                  
end