% Performs all pre-processing steps, such as grid generation, assembly of
% stationary blocks, etc. for the problem solution.

%===============================================================================
%> @file advection/preprocessProblem.m
%>
%> @brief Performs all pre-processing tasks, such as grid generation, assembly
%>        of stationary blocks, etc. for the problem solution.
%===============================================================================
%>
%> @brief Performs all pre-processing steps, such as grid generation, assembly
%>        of stationary blocks, etc. for the problem solution.
%>
%> This routine is called after advection/configureProblem.m.
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
problemData.K    = problemData.g.numT;  % number of triangles
problemData.N    = nchoosek(problemData.p + 2, problemData.p); % number of local DOFs
problemData.Nmu  = problemData.p + 1; % number of local DOFs on Faces
problemData.dt   = problemData.tEnd / problemData.numSteps;

% [K x 3] arrays that mark local edges (E0T) or vertices (V0T) that are
% interior or have a certain boundary type.
problemData.g.markE0Tint  = problemData.generateMarkE0Tint(problemData.g);
problemData.g.markE0TbdrN = problemData.generateMarkE0TbdrN(problemData.g);
problemData.g.markE0TbdrD = problemData.generateMarkE0TbdrD(problemData.g);

% Precompute some repeatedly evaluated fields
problemData.g = computeDerivedGridData(problemData.g);

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.K)

%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct, [problemData.qOrd, problemData.qOrd + 1]);
problemData.basesOnQuad = computeBasesOnQuadEdge(problemData.Nmu, problemData.basesOnQuad, [problemData.qOrd, problemData.qOrd + 1]);

%% Computation of matrices on the reference triangle.
problemData.hatM    = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.hatMmu  = integrateRefEdgeMuMu(problemData.Nmu, problemData.basesOnQuad, problemData.qOrd);
problemData.hatRmu  = integrateRefEdgePhiIntMu([problemData.N, problemData.Nmu], problemData.basesOnQuad, problemData.qOrd);
problemData.hatRphi = integrateRefEdgePhiIntPhiInt(problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.hatG    = integrateRefElemDphiPhiPerQuad(problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.hatS    = integrateRefEdgePhiIntMuPerQuad([problemData.N, problemData.Nmu], problemData.basesOnQuad, problemData.qOrd);

%% Assembly of time-independent global matrices.
problemData.globMphi      = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globRmu       = assembleMatEdgePhiIntMu(problemData.g, problemData.g.markE0Tint, problemData.hatRmu);
problemData.globRphi      = assembleMatEdgePhiIntPhiInt(problemData.g, problemData.g.markE0Tint, problemData.hatRphi);
problemData.globKmuOut    = assembleMatEdgePhiIntMu(problemData.g, problemData.g.markE0TbdrN, problemData.hatRmu)';
problemData.globMmuBar    = assembleMatEdgeMuMu(problemData.g, problemData.g.markE0Tint, problemData.hatMmu);
problemData.globMmuTilde  = assembleMatEdgeMuMu(problemData.g, ~problemData.g.markE0Tint, problemData.hatMmu);
problemData.globP         = problemData.stab .* problemData.globMmuBar + problemData.globMmuTilde;
problemData.globT         = problemData.globRmu';
end % function
