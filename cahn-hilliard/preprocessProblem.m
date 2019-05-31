% Performs all pre-processing steps, such as grid generation, assembly of
% stationary blocks, etc. for the problem solution.

%===============================================================================
%> @file ./cahn-hilliard/preprocessProblem.m
%>
%> @brief Performs all pre-processing tasks, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%===============================================================================
%>
%> @brief Performs all pre-processing steps, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%>
%> This routine is called after cahn-hilliard/configureProblem.m.
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
% %% Store substep handles to avoid calling overhead
% [~, ~, ~, subStepList] = getStepLists();
% problemData.subStepHandles = getStepHandles(problemData.problemName, subStepList);
problemData.stepSuccess = false;
%% Triangulation.
if problemData.Voronoi_tri
  problemData.g = problemData.generateGridData(problemData.hmax, 0, 1, true);
else
  problemData.g = problemData.generateGridData(problemData.hmax);
end
if problemData.isVisGrid,  visualizeGrid(problemData.g);  end
%% Globally constant parameters.
problemData.K           = problemData.g.numT;  % number of triangles
problemData.N           = nchoosek(problemData.p + 2, problemData.p); % number of local DOFs
% problemData.tau         = problemData.tEnd / problemData.numSteps;  % time step size
problemData.t           = 0; % Current time

% [K x 3] arrays that mark local edges (E0T) or vertices (V0T) that are 
% interior or have a certain boundary type.
problemData.g.markE0Tint  = problemData.generateMarkE0Tint(problemData.g);
problemData.g.markE0TbdrD = problemData.generateMarkE0TbdrD(problemData.g);
problemData.g.markE0TbdrN = problemData.generateMarkE0TbdrN(problemData.g);
problemData.g.markV0TbdrD = ismember(problemData.g.V0T, ...
                            problemData.g.V0E(problemData.g.E0T(problemData.g.markE0TbdrD), :));      
%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', problemData.p, problemData.N, problemData.K)
%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(problemData.N, struct);
%% Computation of matrices on the reference triangle.

% Volume contributions
problemData.hatM  = integrateRefElemPhiPhi(problemData.N, problemData.basesOnQuad);
hatL              = integrateRefElemDPhiDPhi(problemData.N, problemData.basesOnQuad);

% Surface contributions
hatSdiag          = integrateRefEdgePhiIntPhiInt(problemData.N, problemData.basesOnQuad);
hatSoffdiag       = integrateRefEdgePhiIntPhiExt(problemData.N, problemData.basesOnQuad);
hatTdiag          = integrateRefEdgeDPhiIntPhiInt(problemData.N, problemData.basesOnQuad);
hatToffdiag       = integrateRefEdgeDPhiIntPhiExt(problemData.N, problemData.basesOnQuad);

% Evaluation at quadrature points
problemData.hatQuadPhiPhi = computeValuesOnQuadElemPhiPhi(problemData.N, problemData.basesOnQuad);
%% Computation of mobility-matrices on the reference triangle.  
hatLmob           = integrateRefElemDPhiDPhiPerQuad(problemData.N, problemData.basesOnQuad);
hatTmobdiag       = integrateRefEdgeDPhiIntPhiIntPerQuad(problemData.N, problemData.basesOnQuad);
hatTmoboffdiag    = integrateRefEdgeDPhiIntPhiExtPerQuad(problemData.N, problemData.basesOnQuad);
hatSmobdiag       = integrateRefEdgePhiIntPhiIntPerQuad(problemData.N, problemData.basesOnQuad);
hatSmoboffdiag    = integrateRefEdgePhiIntPhiExtPerQuad(problemData.N, problemData.basesOnQuad);
problemData.hatLmob = hatLmob;
problemData.hatTmobdiag = hatTmobdiag;
problemData.hatTmoboffdiag = hatTmoboffdiag;
problemData.hatSmobdiag = hatSmobdiag;
problemData.hatSmoboffdiag = hatSmoboffdiag;
%%

% if problemData.isFluxLim || problemData.standardLim
  problemData.basesOnCorners = computeBasesOnCorners(problemData.N);
% end

% if problemData.isFluxLim
  problemData.hatTdiag = hatTdiag;
  problemData.hatToffdiag = hatToffdiag;
  hatBarM = integrateRefElemPhiPhiLumped(problemData.N, problemData.basesOnQuad, problemData.basesOnCorners);
  problemData.hatSdiag = hatSdiag;
  problemData.hatSoffdiag = hatSoffdiag;
% end

if problemData.initialLimiting
  problemData.g.markV0TbdrD = ismember(problemData.g.V0T, ...
                            problemData.g.V0E(problemData.g.E0T(problemData.g.markE0TbdrD), :));
  problemData.basesOnQuad = computeTaylorBasesV0T(problemData.g, problemData.N, problemData.basesOnQuad);
  problemData.globMDiscTaylor = assembleMatElemPhiDiscPhiTaylor(problemData.g, problemData.N, problemData.basesOnQuad);
end

%% Assembly of time-independent global matrices.
% Mass Matrix
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);

% Volume contribution (element diffusion)
globL = assembleMatElemDPhiDPhi(problemData.g, hatL);

% Surface contribution (flux term)
globT  = -1/2 * assembleMatEdgePhiDPhiNu(problemData.g, problemData.g.markE0Tint, ...
                 hatTdiag, hatToffdiag);

% Surface contribution (symmetrization term)
globU  = problemData.eta/2 * assembleMatEdgeDPhiPhiNu(problemData.g, problemData.g.markE0Tint, ...
                 hatTdiag, hatToffdiag);

% Surface contribution (penalty term)
globS  = problemData.sigma * assembleMatEdgePhiPhi(problemData.g, problemData.g.markE0Tint, ...
                 hatSdiag, hatSoffdiag, ones(problemData.K,3));

if problemData.isFluxLim
  problemData.globL = globL;
  problemData.globS = globS;
  problemData.globBarM = assembleMatElemPhiPhi(problemData.g, hatBarM); % Lumped mass matrix
%   problemData.globFe = assembleCellEdgeTermsLimiting(problemData.g, problemData.g.markE0Tint, ...
%       hatTdiag, hatToffdiag, hatSdiag, hatSoffdiag, problemData.eta, problemData.sigma);
end

problemData.globS = globS;
problemData.globA = globL + globT + globU + globS;
end % function