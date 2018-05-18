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
%> This routine is called after swe_2dv/configureProblem.m .
%>
%> This step consists of grid generation, computation of derived
%> data structures, pre-computation of often needed values (e.g., basis
%> functions on quadrature points), and assembly of time-independent matrix
%> blocks.
%>
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields, as provided by swe_2dv/configureProblem.m. 
%>                      @f$[\text{struct}]@f$
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

%% Function handles for time stepping
[~, ~, ~, subStepList] = getStepLists();
problemData.subStepHandles = getStepHandles(problemData.problemName, subStepList);

%% Globally constant parameters.
problemData.N = (problemData.p + 1)^2;  % number of local DOFs on trapezoidals
problemData.barN = problemData.p + 1;  % number of local DOFs on intervals
problemData.tau = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size
problemData.isBottomFriction = ~strcmp(problemData.bottomFriction, 'none');

%% Triangulation.
problemData.g = problemData.generateGrid(problemData.numElem);
problemData.g.g1D = problemData.generateGrid1D(problemData.numElem(1), problemData.g);

%% Additional 2D mesh data: [K x 4] marker for local edges
% Edges that are horizontal or vertical
problemData.g.markE0Th = [ true(problemData.g.numT, 2), false(problemData.g.numT, 2) ];
problemData.g.markE0Tv = ~problemData.g.markE0Th;

% Interior edges
problemData.g.markE0Tint = problemData.generateMarkE0Tint(problemData.g); 

% Vertical boundary edges
problemData.g.markE0Tbdr = sparse(~problemData.g.markE0Tint & problemData.g.markE0Tv); 

% Bottom boundary edges
problemData.g.markE0TbdrBot = sparse(problemData.generateMarkE0TbdrBot(problemData.g));

% Top boundary edges
problemData.g.markE0TbdrTop = sparse(problemData.generateMarkE0TbdrTop(problemData.g));

% Coupling boundary edges
problemData.g.markE0TbdrCoupling = sparse(problemData.generateMarkE0TbdrCoupling(problemData.g));

% Prescribed horizontal velocity
problemData.g.markE0TbdrU = sparse(problemData.generateMarkE0TbdrU(problemData.g));

% Prescribed water height
problemData.g.markE0TbdrH = sparse(problemData.generateMarkE0TbdrH(problemData.g));

% Prescribed diffusion
problemData.g.markE0TbdrQ = sparse(problemData.generateMarkE0TbdrQ(problemData.g));

% Boundary edges with Riemann solver applied
problemData.g.markE0TbdrRiem = sparse(problemData.generateMarkE0TbdrRiem(problemData.g));

assert(~any(problemData.g.markE0Th(:) & problemData.g.markE0TbdrRiem(:)), ...
       'Riemann solver for horizontal boundaries is not implemented!')
assert(nnz(problemData.g.markE0TbdrRiem & ~(problemData.g.markE0TbdrU | problemData.g.markE0TbdrH)) == 0, ...
       'Riemann solver specified for boundaries without Dirichlet data for U or H')

%% Additional 1D mesh data: [barK x 2] marker for local edges
% Interior vertices
problemData.g.g1D.markV0Tint = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0Tint(:, [4 3])) > 0;

% Boundary vertices
problemData.g.g1D.markV0Tbdr = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0Tbdr(:, [4 3])) > 0;

problemData.g.g1D.markV0TbdrU = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0TbdrU(:, [4 3])) > 0;
problemData.g.g1D.markV0TbdrH = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0TbdrH(:, [4 3])) > 0;

% Boundary vertices with Riemann solver applied
problemData.g.g1D.markV0TbdrRiem = problemData.g.g1D.markT2DT.' * double(problemData.g.markE0TbdrRiem(:, [4 3])) > 0;

%% Function handles for problem-specific functions
problemData.fn_adaptFreeSurface = getFunctionHandle([problemData.problemName filesep 'adaptFreeSurface']);
problemData.fn_assembleMatEdgeTetraHorizPhiPhiNuBottomUp = getFunctionHandle([problemData.problemName filesep 'assembleMatEdgeTetraHorizPhiPhiNuBottomUp']);
problemData.fn_assembleMatEdgeTetraVertPhiPhiFuncDisc1DNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatEdgeTetraVertPhiPhiFuncDisc1DNuHeight']);
problemData.fn_assembleMatEdgeTetraVertPhiIntPhiIntFuncDisc1DIntNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatEdgeTetraVertPhiIntPhiIntFuncDisc1DIntNuHeight']);
problemData.fn_assembleMatElem1DDphiPhiFuncDiscHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatElem1DDphiPhiFuncDiscHeight']);
problemData.fn_assembleMatV0T1DPhiPhiFuncDiscNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatV0T1DPhiPhiFuncDiscNuHeight']);
problemData.fn_assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight']);
problemData.fn_assembleVecEdgeTetraVertPhiIntFuncContHeightNu = getFunctionHandle([problemData.problemName filesep 'assembleVecEdgeTetraVertPhiIntFuncContHeightNu']);
problemData.fn_assembleVecEdgeTetraVertPhiIntFuncDiscIntHeightNu = getFunctionHandle([problemData.problemName filesep 'assembleVecEdgeTetraVertPhiIntFuncDiscIntHeightNu']);
problemData.fn_assembleVecV0T1DPhiIntFuncContNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleVecV0T1DPhiIntFuncContNuHeight']);
problemData.fn_assembleVecV0T1DPhiIntFuncDiscIntNuHeight = getFunctionHandle([problemData.problemName filesep 'assembleVecV0T1DPhiIntFuncDiscIntNuHeight']);

%% Configuration output.
fprintf('-------------------------------------------------------------------------------------------\n');
fprintf('Running problem "%s" with testcase "%s".\n', problemData.problemName, problemData.testcase);
fprintf('Computing with polynomial order %d (%d local DOFs) on %d x %d (%d) trapezoids.\n', ...
        problemData.p, problemData.N, problemData.numElem(1), problemData.numElem(2), problemData.g.numT);
fprintf('%d time steps from t = %g to %g.\n', problemData.numSteps, problemData.t0, problemData.tEnd);
fprintf('-------------------------------------------------------------------------------------------\n');

%% Lookup tables for basis function.
assert(problemData.qOrd <= problemData.qOrdMax, 'Maximum requested order of quadrature rule must be greater than or equal to qOrd')
problemData.basesOnQuad1D = computeBasesOnQuad1D(problemData.p, struct, problemData.qOrd : problemData.qOrdMax+1);
problemData.basesOnQuad2D = computeBasesOnQuadTensorProduct(problemData.p, struct, problemData.qOrd : problemData.qOrdMax+1);

%% Computation of matrices on the reference element.
problemData.hatM = integrateRefElemTetraPhiPhi(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatG = integrateRefElemTetraDphiPhiPhi(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatH = integrateRefElemTetraDphiPhi(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatQdiag = integrateRefEdgePhiIntPhiInt(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatQoffdiag = integrateRefEdgePhiIntPhiExt(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatRdiag = integrateRefEdgePhiIntPhiIntPhiInt(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatRoffdiag = integrateRefEdgePhiIntPhiExtPhiExt(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
problemData.hatSdiag = integrateRefEdgePhiIntPerQuad(problemData.N, problemData.basesOnQuad2D, problemData.qOrd);

problemData.hatBarM = integrateRefElem1DPhiPhi(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.hatBarG = integrateRefElem1DDphiPhiPhiPerQuad(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.hatBarSdiag = integrateRefEdgeTetraPhi1DIntPerQuad(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.hatBarPdiag = computePhiIntPhiIntPhiIntV0T1D(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.hatBarPoffdiag = computePhiIntPhiExtPhiExtV0T1D(problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);

problemData.hatVeeH = integrateRefElemTetraDphiPhi1D([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.hatVeeQdiag = integrateRefEdgeTetraPhiIntPhi1DInt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.hatVeeQoffdiag = integrateRefEdgeTetraPhiIntPhi1DExt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.hatVeePdiag = integrateRefEdgeTetraPhiIntPhiIntPhi1DInt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);
problemData.hatVeePoffdiag = integrateRefEdgeTetraPhiIntPhiExtPhi1DExt([problemData.N problemData.barN], problemData.qOrd, problemData.basesOnQuad2D, problemData.basesOnQuad1D);

%% One-dimensional mass matrix in free-surface equation 
problemData.globBarM = assembleMatElemPhiPhi(problemData.g.g1D, problemData.hatBarM);

%% Empty vectors and matrices for coupled problem
problemData.globJuCoupling = { sparse(problemData.g.numT * problemData.N, 1), sparse(problemData.g.numT * problemData.N, 1) };
problemData.globJwCoupling = sparse(problemData.g.numT * problemData.N, 1);
problemData.globJuuCoupling = { sparse(problemData.g.numT * problemData.N, 1), sparse(problemData.g.numT * problemData.N, 1) };
end % function