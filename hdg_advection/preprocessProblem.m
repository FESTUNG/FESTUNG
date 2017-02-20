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

problemData.cDiscRK = cell( problemData.tabRK.s, 1);

%% Globally constant parameters.
problemData.K           = problemData.g.numT;  % number of triangles
problemData.N           = nchoosek(problemData.p + 2, problemData.p); % number of local DOFs
problemData.Nlambda     = problemData.p + 1; % number of local DOFs on Faces
problemData.tau         = problemData.tEnd / problemData.numSteps;  % time step size

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
problemData.g.flipArray = generateFlipArray( problemData.g );
% Choose a block size for the local solves if we want 'true' local solves
if ( problemData.isTrueLocalSolve  == true )
    K = problemData.K;
    problemData.localSolveBlockSize = min(K, 16);
    if (mod(K, problemData.localSolveBlockSize) ~= 0)
        problemData.localSolveBlockSize = 1;
        warning('Block size does not fit the problem. Reset block size to 1');
    end
    assert( mod(K, problemData.localSolveBlockSize) == 0, ...
        'Block size does not fit to problem size!');
end

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.K)
%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct);
problemData.basesOnGamma = computeBasesOnGamma(problemData.Nlambda, struct);

%% Computation of matrices on the reference triangle.
problemData.hatM              = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad);

%% Computation of HDG matrices on the reference triangle.
%Hybrid mass matrix
problemData.hatMlambda = integrateRefEdgeMuMu(problemData.Nlambda, problemData.basesOnGamma);

% Integrals on edges
problemData.hatRlambda = integrateRefEdgePhiIntMu(problemData.N, problemData.Nlambda, problemData.basesOnQuad, problemData.basesOnGamma);
problemData.hatRphi    = integrateRefEdgePhiIntPhiInt(problemData.N, problemData.basesOnQuad);
% Precomputations for term II
problemData.hatGbarOnQuad = integrateRefElemDphiPhiFlux(problemData.N, problemData.basesOnQuad);
% Precomputations for III.1
problemData.hatSbarOnQuad = integrateRefEdgeMuPhiIntFlux(problemData.N, problemData.Nlambda, problemData.basesOnQuad, problemData.basesOnGamma);

%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);

%Term III.2 and Term III.5
problemData.globRD  = assembleMatEdgeMuPhiInt( problemData.g, problemData.g.markE0TbdrN, problemData.hatRlambda );

%Term III.2, we use the assembly routine above
%Interior
problemData.globRlambdaBar = assembleMatEdgeMuPhiInt( problemData.g, problemData.g.markE0Tint, problemData.hatRlambda );
%Term III.5, we use the assembly routine above
% Exterior
problemData.globRlambdaHat = assembleMatEdgeMuPhiInt( problemData.g, ~problemData.g.markE0Tint, problemData.hatRlambda );
% Assemble total matrix
problemData.globRlambda = problemData.globRlambdaBar + problemData.globRlambdaHat;

%Term III.3 WIP
problemData.globRphi = assembleMatEdgePhiIntPhiIntHybrid( problemData.g, problemData.hatRphi );
%Term VI.2
problemData.globRgamma = assembleMatEdgeMuPhiInt( problemData.g, problemData.g.markE0TbdrN, problemData.hatRlambda );
problemData.globRgamma = problemData.globRgamma';

%Term V.1
problemData.globMint  = assembleMatEdgeMuMu(problemData.g,  problemData.g.markE0Tint , problemData.hatMlambda);
%Term VI.1
problemData.globMext  = assembleMatEdgeMuMu(problemData.g, ~problemData.g.markE0Tint , problemData.hatMlambda);
%Matrix P
problemData.globP = problemData.stab .* problemData.globMint + problemData.globMext ;

%Term V.2 g, markE0Tbdr, refEdgePhiIntMu
problemData.globU = assembleMatEdgeMuPhiInt( problemData.g, problemData.g.markE0Tint, problemData.hatRlambda );
problemData.globU = problemData.globU';

if ( problemData.showWaitBar == true )
    problemData.waitBar = waitbar( 0, 'Simulation progress');
end
end % function
