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
%> This routine is called after darcy_2dv/configureProblem.m.
%>
%> This step consists of grid generation, computation of derived
%> data structures, pre-computation of often needed values (e.g., basis
%> functions on quadrature points), and assembly of time-independent matrix
%> blocks.
%>
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields, as provided by darcy_2dv/configureProblem.m. 
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

%% Triangulation.
problemData.g = problemData.generateGrid(problemData.numElem);
if problemData.isVisGrid, visualizeGridQuadri(problemData.g); end

%% Globally constant parameters.
problemData.N = (problemData.p + 1)^2;  % number of local DOFs
problemData.tau = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size

%% Additional mesh data
% [K x 3] mark local edges that are interior or boundary
problemData.g.markE0Tint = problemData.generateMarkE0Tint(problemData.g); 
problemData.g.markE0TbdrN = sparse(problemData.generateMarkE0TbdrN(problemData.g));
problemData.g.markE0TbdrD = sparse(problemData.generateMarkE0TbdrD(problemData.g));
problemData.g.markE0TbdrCoupling = sparse(problemData.generateMarkE0TbdrCoupling(problemData.g));

%% Configuration output.
fprintf('-------------------------------------------------------------------------------------------\n');
fprintf('Running testcase "%s".\n', problemData.testcase);
fprintf('Computing with polynomial order %d (%d local DOFs) on %d x %d (%d) trapezoids.\n', ...
        problemData.p, problemData.N, problemData.numElem(1), problemData.numElem(2), problemData.g.numT);
if problemData.isStationary
  fprintf('Computing stationary solution.\n');
else
  fprintf('%d time steps from t = %g to %g.\n', problemData.numSteps, problemData.t0, problemData.tEnd);
end % if
fprintf('-------------------------------------------------------------------------------------------\n');

%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuadTensorProduct(problemData.p, struct, problemData.qOrd : problemData.qOrdMax+1);

%% Computation of matrices on the reference element.
problemData.hatM = integrateRefElemQuadriPhiPhi(problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.hatG = integrateRefElemQuadriDphiPhiPhi(problemData.N, problemData.basesOnQuad, problemData.qOrd);
hatH = integrateRefElemQuadriDphiPhi(problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.hatRdiag = integrateRefEdgePhiIntPhiIntPhiInt(problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.hatRoffdiag = integrateRefEdgePhiIntPhiExtPhiExt(problemData.N, problemData.basesOnQuad, problemData.qOrd);
hatSdiag = integrateRefEdgePhiIntPhiInt(problemData.N, problemData.basesOnQuad, problemData.qOrd);
hatSoffdiag = integrateRefEdgePhiIntPhiExt(problemData.N, problemData.basesOnQuad, problemData.qOrd);

%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
problemData.globH = assembleMatElemDphiPhi(problemData.g, hatH);
problemData.globQ = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint, hatSdiag, hatSoffdiag);
problemData.globQN = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrN, hatSdiag);
problemData.globS = problemData.eta * assembleMatEdgePhiPhi(problemData.g, problemData.g.markE0Tint, hatSdiag, hatSoffdiag, ones(problemData.g.numT, 4));
if problemData.isJumpCoupling
  problemData.globSD = problemData.eta * assembleMatEdgePhiIntPhiInt(problemData.g, problemData.g.markE0TbdrD | problemData.g.markE0TbdrCoupling, hatSdiag, ones(problemData.g.numT, 4));
else
  problemData.globSD = problemData.eta * assembleMatEdgePhiIntPhiInt(problemData.g, problemData.g.markE0TbdrD, hatSdiag, ones(problemData.g.numT, 4));
end % if

if ~problemData.isStationary
  problemData.sysW = [ sparse(2 * problemData.g.numT * problemData.N, 3 * problemData.g.numT * problemData.N) ; ...
                       sparse(problemData.g.numT * problemData.N, 2 * problemData.g.numT * problemData.N), problemData.globM ];
end % if

%% Empty vectors for coupled problem
problemData.globJcouple = { sparse(problemData.g.numT * problemData.N, 1), sparse(problemData.g.numT * problemData.N, 1) };
problemData.globKcouple = sparse(problemData.g.numT * problemData.N, 1);
end % function