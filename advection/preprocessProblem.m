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
%% Store substep handles to avoid calling overhead
[~, ~, ~, subStepList] = getStepLists();
problemData.subStepHandles = getStepHandles(problemData.problemName, subStepList);
%% Triangulation.
problemData.g = problemData.generateGridData(problemData.hmax);
if problemData.isVisGrid
  if problemData.isQuadri
    visualizeGridQuadri(problemData.g);
  else
    visualizeGrid(problemData.g);  
  end
end
%% Globally constant parameters.
problemData.K = problemData.g.numT;  % number of triangles
if problemData.isQuadri
  problemData.N = (problemData.p + 1) * (problemData.p + 1);
else
  problemData.N = nchoosek(problemData.p + 2, problemData.p); % number of local DOFs
end
problemData.tau = problemData.tEnd / problemData.numSteps;  % time step size

% [K x 3 or 4] arrays that mark local edges (E0T) or vertices (V0T) that are 
% interior or have a certain boundary type.
problemData.g.markE0Tint  = problemData.generateMarkE0Tint(problemData.g);
problemData.g.markE0TbdrN = problemData.generateMarkE0TbdrN(problemData.g);
problemData.g.markE0TbdrD = problemData.generateMarkE0TbdrD(problemData.g);
problemData.g.markV0TbdrD = ismember(problemData.g.V0T, ...
                            problemData.g.V0E(problemData.g.E0T(problemData.g.markE0TbdrD), :)); 

% Precompute some repeatedly evaluated fields                          
if ~problemData.isQuadri
  problemData.g = computeDerivedGridData(problemData.g);
end
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.K)
%% Lookup table for basis function.
if problemData.isQuadri
  problemData.basesOnQuad = computeBasesOnQuadTensorProduct(problemData.p, struct, problemData.qOrd : problemData.qOrdMax+1);
else
  problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct);
end
if problemData.isSlopeLim
  problemData.basesOnQuad = computeTaylorBasesV0T(problemData.g, problemData.N, problemData.basesOnQuad);
end % if
%% Computation of matrices on the reference triangle.
if problemData.isQuadri
  problemData.hatM = integrateRefElemQuadriPhiPhi(problemData.N, problemData.basesOnQuad, problemData.qOrd);
  problemData.hatG = integrateRefElemQuadriDphiPhiPhi(problemData.N, problemData.basesOnQuad, problemData.qOrd);
  problemData.hatRdiagOnQuad = integrateRefEdgeQuadriPhiIntPhiIntPerQuad(problemData.N, problemData.basesOnQuad, problemData.qOrd);
  problemData.hatRoffdiagOnQuad = integrateRefEdgeQuadriPhiIntPhiExtPerQuad(problemData.N, problemData.basesOnQuad, problemData.qOrd);
else
  problemData.hatM              = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad);
  problemData.hatG              = integrateRefElemDphiPhiPhi(problemData.N, problemData.basesOnQuad);
  problemData.hatRdiagOnQuad    = integrateRefEdgePhiIntPhiIntPerQuad(problemData.N, problemData.basesOnQuad);
  problemData.hatRoffdiagOnQuad = integrateRefEdgePhiIntPhiExtPerQuad(problemData.N, problemData.basesOnQuad);
end
%% Assembly of time-independent global matrices.
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);
if problemData.isSlopeLim
  globMTaylor = assembleMatElemPhiTaylorPhiTaylor(problemData.g, problemData.N);
  problemData.globMDiscTaylor = assembleMatElemPhiDiscPhiTaylor(problemData.g, problemData.N, problemData.basesOnQuad);
  problemData.globMCorr = spdiags(1./diag(globMTaylor), 0, problemData.K * problemData.N, problemData.K * problemData.N) * globMTaylor;
end % if
end % function