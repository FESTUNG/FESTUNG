% Performs all pre-processing steps, such as grid generation, assembly of
% stationary blocks, etc. for the problem solution.

%===============================================================================
%> @file template/preprocessProblem.m
%>
%> @brief Performs all pre-processing tasks, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%===============================================================================
%>
%> @brief Performs all pre-processing steps, such as grid generation, assembly 
%>        of stationary blocks, etc. for the problem solution.
%>
%> This routine is called after template/configureProblem.m.
%>
%> Typically, this step consists of grid generation, computation of derived
%> data structures, pre-computing often needed values (e.g., basis
%> functions on quadrature points), or assembly of time-independent matrix
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
p = problemData.p;

%% Output preprocessing
if ~isdir('output')
  mkdir('output');
end % if

%% Triangulation
X1 = [0 1 1 0]; X2 = [0 0 1 1];
problemData.g = domainHierarchy(X1, X2, problemData.hmax, problemData.refinement);
problemData.g.idE = zeros(problemData.g.numE,1);
problemData.g.idE(problemData.g.baryE(:, 2) == 0) = 1; % south
problemData.g.idE(problemData.g.baryE(:, 1) == 1) = 4; % east
problemData.g.idE(problemData.g.baryE(:, 2) == 1) = 1; % north
problemData.g.idE(problemData.g.baryE(:, 1) == 0) = 1; % west
problemData.g.idE0T = problemData.g.idE(problemData.g.E0T);

%% Globally constant parameters
problemData.K = problemData.g.numT; % number of triangles
problemData.N = nchoosek(p + 2, p); % number of local DOFs
K = problemData.K;
N = problemData.N;

problemData.g.markE0Tint = problemData.g.idE0T == 0; % [K x 3] mark local edges that are interior
problemData.g.markE0TbdrL = problemData.g.idE0T == 1; % [K x 3] mark local edges on the open sea boundary
problemData.g.markE0TbdrRA = problemData.g.idE0T == 2; % [K x 3] mark local edges on the open sea boundary
problemData.g.markE0TbdrRI = problemData.g.idE0T == 3; % [K x 3] mark local edges on the open sea boundary
problemData.g.markE0TbdrOS = problemData.g.idE0T == 4; % [K x 3] mark local edges on the open sea boundary
riverBdrs = ~isequal(problemData.g.markE0TbdrRI, zeros(K, 3));
problemData.g = computeDerivedGridData(problemData.g);

problemData.outputFrequency = max(floor(problemData.numSteps / problemData.outputCount), 1);

%% Configuration output.
fprintf('Computing with polynomial order %d (%d local DOFs) on %d triangles.\n', p, N, K);

%% Lookup table for basis function.
problemData.basesOnQuad = computeBasesOnQuad(N, struct);

%% System matrix for correction of min value exceedence.
problemData.sysMinValueCorrection = [ phi(1,0,0) phi(1,1,0) phi(1,0,1) ; ...
                                      phi(2,0,0) phi(2,1,0) phi(2,0,1) ; ...
                                      phi(3,0,0) phi(3,1,0) phi(3,0,1) ];

%% Computation of matrices on the reference triangle.
problemData.refElemPhiPhiPhi = integrateRefElemPhiPhiPhi(N, problemData.basesOnQuad);
problemData.refElemDphiLinPhiPhi = integrateRefElemDphiLinPhiPhi(N, problemData.basesOnQuad);
problemData.refElemDphiPhi = integrateRefElemDphiPhi(N, problemData.basesOnQuad);
problemData.refElemPhiPhi = integrateRefElemPhiPhi(N, problemData.basesOnQuad);
problemData.refEdgePhiIntPhiInt = integrateRefEdgePhiIntPhiInt(N, problemData.basesOnQuad);
problemData.refEdgePhiIntPhiExt = integrateRefEdgePhiIntPhiExt(N, problemData.basesOnQuad);
problemData.refElemDphiPerQuad = integrateRefElemDphiPerQuad(N, problemData.basesOnQuad);
problemData.refEdgePhiIntPerQuad = integrateRefEdgePhiIntPerQuad(N, problemData.basesOnQuad);
problemData.refEdgePhiIntPhiIntPerQuad = integrateRefEdgePhiIntPhiIntPerQuad(N, problemData.basesOnQuad);

%% L2 projections of time-independent algebraic coefficients.
problemData.zbDisc = projectFuncCont2DataDisc(problemData.g, problemData.zbCont, 2, eye(3), computeBasesOnQuad(3, struct)); 
fcDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fcCont(x1,x2), 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);

% Visualization of coefficients
if any(ismember(problemData.outputList, 'zb'))
  dataLagr = projectDataDisc2DataLagr(problemData.zbDisc);
  visualizeDataLagr(problemData.g, dataLagr, 'z_b', ['output/' problemData.name '_zb'], 0, problemData.outputTypes);
end % if
if any(ismember(problemData.outputList, 'fc'))
  dataLagr = projectDataDisc2DataLagr(fcDisc);
  visualizeDataLagr(problemData.g, dataLagr, 'f_c', ['output/' problemData.name '_fc'], 0, problemData.outputTypes);
end % if

%% Assembly of time-independent global matrices corresponding to linear contributions.
globD = assembleMatElemPhiPhiFuncDisc(problemData.g, problemData.refElemPhiPhiPhi, fcDisc);
globG = assembleMatElemPhiPhiFuncDiscLin(problemData.g, problemData.refElemDphiLinPhiPhi, problemData.zbDisc);
globH = assembleMatElemDphiPhi(problemData.g, problemData.refElemDphiPhi);
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.refElemPhiPhi);
globQ = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint, problemData.refEdgePhiIntPhiInt, problemData.refEdgePhiIntPhiExt, problemData.g.areaNuE0Tint);
globQOS = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrOS, problemData.refEdgePhiIntPhiInt, problemData.g.areaNuE0TbdrOS);
globQRA = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrRA, problemData.refEdgePhiIntPhiInt, problemData.g.areaNuE0TbdrRA);

problemData.sysW = blkdiag(problemData.globM, problemData.globM, problemData.globM);
problemData.linearTerms = [               sparse(K*N,K*N), globQ{1} + globQOS{1} + globQRA{1} - globH{1},  globQ{2} + globQOS{2} + globQRA{2} - globH{2}; ...
                            problemData.gConst * globG{1},                               sparse(K*N,K*N),                                          -globD; ...
                            problemData.gConst * globG{2},                                         globD,                                 sparse(K*N,K*N) ];
%% Assembly of time-independent global matrices corresponding to non-linear contributions.

end % function