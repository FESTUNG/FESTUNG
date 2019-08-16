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
%> This routine is called after adv-diff-reac/configureProblem.m.
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
%> @copyright 2014-2019 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
problemData.N = nchoosek(problemData.p + 2, problemData.p);  % number of local DOFs
problemData.dt = (problemData.tEnd - problemData.t0) / problemData.numSteps;  % time step size
%% Additional mesh data
% [K x 3] arrays that mark local edges that are interior or have a certain boundary type.
problemData.g.markE0Tint = problemData.generateMarkE0Tint(problemData.g);
problemData.g.markE0TbdrN = problemData.generateMarkE0TbdrN(problemData.g);
problemData.g.markE0TbdrD = problemData.generateMarkE0TbdrD(problemData.g);
% problemData.g = computeDerivedGridData(problemData.g);
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', ...
        problemData.p, problemData.N, problemData.g.numT)
if problemData.isIP
  fprintf('Using interior penalty dG discretization with theta=%.2f and eta=%.2f.\n', ...
          problemData.symparam, problemData.penparam)
else
  fprintf('Using LDG discretization with eta=%f.\n', problemData.penparam)
end % if
%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct);
%% Computation of matrices on the reference triangle.
problemData.refElemPhiPhi = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad);
problemData.refElemPhiPhiPerQuad = integrateRefElemPhiPhiPerQuad(problemData.N, problemData.basesOnQuad);
problemData.refElemDphiPhiPerQuad = integrateRefElemDphiPhiPerQuad(problemData.N, problemData.basesOnQuad);
refEdgePhiIntPhiInt = integrateRefEdgePhiIntPhiInt(problemData.N, problemData.basesOnQuad);
refEdgePhiIntPhiExt = integrateRefEdgePhiIntPhiExt(problemData.N, problemData.basesOnQuad);
problemData.refEdgePhiIntPhiIntPerQuad = integrateRefEdgePhiIntPhiIntPerQuad(problemData.N, problemData.basesOnQuad);
problemData.refEdgePhiIntPhiExtPerQuad = integrateRefEdgePhiIntPhiExtPerQuad(problemData.N, problemData.basesOnQuad);
%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.refElemPhiPhi);
if problemData.penparam == 0
  problemData.globBjmp = sparse(problemData.g.numT * problemData.N, problemData.g.numT * problemData.N);
  problemData.globJmp = zeros(problemData.g.numT * problemData.N, 1);
else
  problemData.globBjmp = assembleMatEdgePhiPhi(problemData.g, ~problemData.g.markE0TbdrN, ...
    refEdgePhiIntPhiInt, refEdgePhiIntPhiExt, ones(problemData.g.numT, 3));
end % if
%% Computation of reference blocks and global matrices specific to IP-dG or LDG
if problemData.isIP  % use IP discretization of diffusion
  problemData.refElemDphiDphiPerQuad = integrateRefElemDphiDphiPerQuad(problemData.N, problemData.basesOnQuad);
  problemData.refEdgeDphiIntPhiIntPerQuad = integrateRefEdgeDphiIntPhiIntPerQuad(problemData.N, problemData.basesOnQuad);
  problemData.refEdgeDphiIntPhiExtPerQuad = integrateRefEdgeDphiIntPhiExtPerQuad(problemData.N, problemData.basesOnQuad);
  if problemData.symparam == 0
    problemData.globBsym = sparse(problemData.g.numT * problemData.N, problemData.g.numT * problemData.N);
    problemData.globJsym = zeros(problemData.g.numT * problemData.N, 1);
  end % if
else                 % use LDG discretization of diffusion
  refElemDphiPhi = integrateRefElemDphiPhi(problemData.N, problemData.basesOnQuad);
  problemData.globAq = assembleMatElemDphiPhi(problemData.g, refElemDphiPhi);
  problemData.globBq = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint, ...
                                               refEdgePhiIntPhiInt, refEdgePhiIntPhiExt);
  problemData.globBqN = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrN, refEdgePhiIntPhiInt);
end
end