% Performs all pre-processing steps, such as grid generation, assembly of
% stationary blocks, etc. for the problem solution.

%===============================================================================
%> @file
%>
%> @brief Performs all pre-processing tasks, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%===============================================================================
%>
%> @brief Performs all pre-processing steps, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%>
%> This routine is called after diffusion/configureProblem.m.
%>
%> This step consists of grid generation, computation of derived
%> data structures, pre-computation of often needed values (e.g., basis
%> functions on quadrature points), and assembly of time-independent matrix
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
function problemData = preprocessProblem(problemData)
%% Triangulation.
problemData.g = problemData.generateGridData(problemData.hmax);
if problemData.isVisGrid,  visualizeGrid(problemData.g);  end
%% Globally constant parameters.
problemData.K = problemData.g.numT;  % number of triangles
problemData.N = nchoosek(problemData.p + 2, problemData.p);  % number of local DOFs
problemData.tau = problemData.tEnd / problemData. numSteps;  % time step size
%% Additional mesh data
% [K x 3] arrays that mark local edges that are interior or have a certain boundary type.
problemData.g.markE0Tint  = problemData.generateMarkE0Tint(problemData.g);
problemData.g.markE0TbdrN = problemData.generateMarkE0TbdrN(problemData.g);
problemData.g.markE0TbdrD = problemData.generateMarkE0TbdrD(problemData.g);
problemData.g = computeDerivedGridData(problemData.g);
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.K)
%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct);
%% Computation of matrices on the reference triangle.
problemData.hatM = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad);
problemData.hatG = integrateRefElemDphiPhiPhi(problemData.N, problemData.basesOnQuad);
problemData.hatH = integrateRefElemDphiPhi(problemData.N, problemData.basesOnQuad);
problemData.hatRdiag = integrateRefEdgePhiIntPhiIntPhiInt(problemData.N, problemData.basesOnQuad);
problemData.hatRoffdiag = integrateRefEdgePhiIntPhiExtPhiExt(problemData.N, problemData.basesOnQuad);
problemData.hatSdiag = integrateRefEdgePhiIntPhiInt(problemData.N, problemData.basesOnQuad);
problemData.hatSoffdiag = integrateRefEdgePhiIntPhiExt(problemData.N, problemData.basesOnQuad);
%% Assembly of time-independent global matrices.
problemData.globM  = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globH  = assembleMatElemDphiPhi(problemData.g, problemData.hatH);
problemData.globQ  = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatSdiag, problemData.hatSoffdiag);
problemData.globQN = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrN, ...
                      problemData.hatSdiag, problemData.g.areaNuE0TbdrN);
problemData.globS  = problemData.eta * assembleMatEdgePhiPhi(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatSdiag, problemData.hatSoffdiag);
problemData.globSD = problemData.eta * assembleMatEdgePhiIntPhiInt(problemData.g, ...
                      problemData.g.markE0TbdrD, problemData.hatSdiag, ones(problemData.K, 3));
if ~problemData.isStationary
  problemData.sysW = [ sparse(2 * problemData.K * problemData.N, 3 * problemData.K * problemData.N) ; ...
                       sparse(problemData.K * problemData.N, 2 * problemData.K * problemData.N), problemData.globM ];
end % if
end % function