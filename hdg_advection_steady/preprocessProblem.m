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
% problemData.tau  = problemData.tEnd / problemData.numSteps;  % time step size

% [K x 3] arrays that mark local edges (E0T) or vertices (V0T) that are
% interior or have a certain boundary type.
problemData.g.markE0Tint  = problemData.generateMarkE0Tint(problemData.g);
problemData.g.markE0TbdrN = problemData.generateMarkE0TbdrN(problemData.g);
problemData.g.markE0TbdrD = problemData.generateMarkE0TbdrD(problemData.g);
problemData.g.markV0TbdrD = ismember(problemData.g.V0T, ...
    problemData.g.V0E(problemData.g.E0T(problemData.g.markE0TbdrD), :));

% Precompute some repeatedly evaluated fields
problemData.g = computeDerivedGridData(problemData.g);

%% HDG related configuration
problemData.g.markSideE0T = generateMarkSideE0T( problemData.g );
% Choose a block size for the local solves if we want 'true' local solves
if ( problemData.isTrueLocalSolve  == true )
    problemData.localSolveBlockSize = determineLocalSolveBlockSize( problemData.K );
end

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.K)
%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct);
problemData.basesOnGamma = computeBasesOnGamma(problemData.Nmu, struct);

%% Computation of matrices on the reference triangle.
problemData.hatM              = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad);

%% Computation of HDG matrices on the reference triangle/edge2.
%Hybrid mass matrix
problemData.hatMmu = integrateRefEdgeMuMu(problemData.Nmu, problemData.basesOnGamma);

% Integrals on edges
problemData.hatRmu  = integrateRefEdgePhiIntMu(problemData.N, problemData.Nmu, problemData.basesOnQuad, problemData.basesOnGamma);
problemData.hatRphi = integrateRefEdgePhiIntPhiInt(problemData.N, problemData.basesOnQuad);
% Precomputations for term II
problemData.hatG = integrateRefElemDphiPhiFlux(problemData.N, problemData.basesOnQuad);
% Precomputations for III.1
problemData.hatS = integrateRefEdgeMuPhiIntFlux(problemData.N, problemData.Nmu, problemData.basesOnQuad, problemData.basesOnGamma);

%% Assembly of time-independent global matrices.
problemData.globMphi = assembleMatElemPhiPhi(problemData.g, problemData.hatM);

%We use the assembly routine above
%Interior
problemData.globRmu = assembleMatEdgeMuPhiInt( problemData.g, problemData.g.markE0Tint, problemData.hatRmu );
%
problemData.globRphi = assembleMatEdgePhiIntPhiIntHybrid( problemData.g, problemData.hatRphi );
%
problemData.globKmuOut = assembleMatEdgeMuPhiInt( problemData.g, problemData.g.markE0TbdrN, problemData.hatRmu );
problemData.globKmuOut = problemData.globKmuOut';

%
problemData.globMmuBar  = assembleMatEdgeMuMu(problemData.g,  problemData.g.markE0Tint , problemData.hatMmu);
%
problemData.globMmuTilde  = assembleMatEdgeMuMu(problemData.g, ~problemData.g.markE0Tint , problemData.hatMmu);
%
problemData.globP = problemData.stab .* problemData.globMmuBar + problemData.globMmuTilde;

%
problemData.globT = assembleMatEdgeMuPhiInt( problemData.g, problemData.g.markE0Tint, problemData.hatRmu );
problemData.globT = problemData.globT';

if ( problemData.showWaitBar == true )
    problemData.waitBar = waitbar( 0, 'Simulation progress');
end
end % function
