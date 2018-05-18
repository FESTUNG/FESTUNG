% Adapts the mesh to match the free surface elevation and reassembles static matrices.

%===============================================================================
%> @file
%>
%> @brief Adapts the mesh to match the free surface elevation and reassembles
%>        static matrices.
%===============================================================================
%>
%> @brief Adapts the mesh to match the free surface elevation and reassembles
%>        static matrices.
%>
%> This routine evaluates the water height @f$h@f$ in each surface vertex of the
%> mesh and derives the smoothed mesh height @f$H_\mathrm{s}@f$ from it
%> (see @ref RRAFK2018).
%> In the case of a difference between the current height of the mesh and the
%> the computed smoothed mesh height it adapts the @f$x^2@f$ corrdinates of the 
%> surface vertices, updates mesh data structures accordingly, and calls 
%> assembleStaticMatrices() to trigger re-assembly of dependent matrices.
%>
%> @param  problemData  A struct with problem parameters, basis functions and
%>                      representation matrices of the current solution
%>                      @f$[\text{struct}]@f$
%> @param  force        Force re-assembly of static matrices 
%>                      @f$[1\times1\text{ bool}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
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
function [problemData] = adaptFreeSurface(problemData, force)
if nargin < 2
  force = false;
end % if

% cDisc{1} in vertices of surface mesh
problemData.hV0T1D = problemData.cDiscRK{end, 1} * problemData.basesOnQuad1D.phi0D{problemData.qOrd}; 

% smoothed height in vertices of surface mesh
hSmoothV0T1D = problemData.hV0T1D;
for n = 1 : 2
  hSmoothV0T1D(:, n) = (1 ./ (1 + sum(problemData.g.g1D.markV0TV0T{n}, 2))) .* ...
                      (hSmoothV0T1D(:, n) + problemData.g.g1D.markV0TV0T{n} * problemData.hV0T1D(:, 3-n));
end % for n

% surface elevation in vertices of surface mesh
xiSmoothV0T1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), [1 2], 2) + hSmoothV0T1D;

% Check for negative or complex height
if ~isreal(hSmoothV0T1D) || any(hSmoothV0T1D(:) < 0)
  error('Negative or complex height');
end % if

% check for necessity of adaptation
markModifiedV0T1D = problemData.g.g1D.coordV0T(:, :, 2) ~= xiSmoothV0T1D;
if force || any(markModifiedV0T1D(:))
  % Change 1D coordinates
  idxV1D = problemData.g.g1D.V0T(markModifiedV0T1D);
  problemData.g.g1D.coordV(idxV1D,2) = xiSmoothV0T1D(markModifiedV0T1D);
  for n = 1 : 2
    problemData.g.g1D.coordV0T(markModifiedV0T1D(:, n), n, 2) = xiSmoothV0T1D(markModifiedV0T1D(:, n), n);
  end % for n
  
  % Change 2D coordinates
  idxV2D = problemData.g.g1D.idxV2D0V(idxV1D);
  problemData.g.coordV(idxV2D, 2) = xiSmoothV0T1D(markModifiedV0T1D);
  for k = 3 : 4
    problemData.g.coordV0T(:, k, 2) = problemData.g.coordV(problemData.g.V0T(:, k), 2);
  end % for k
  
  % Re-generate coordinate-dependent grid data
  problemData.g = problemData.g.generateCoordDependGridData(problemData.g);
  
  % Determine need for projection
  isProjectData = all(isfield(problemData, {'globM', 'cDiscRK'}));
  if isProjectData, oldGlobM = problemData.globM; end
  
  % Re-assemble static matrices
  problemData = assembleStaticMatrices(problemData);
  
  % Project velocity to new grid
  if isProjectData
    cSys = [ reshape(problemData.cDiscRK{end, 2}.', [], 1), reshape(problemData.cDiscRK{end, 3}.', [], 1) ];
    cSys = problemData.globM \ (oldGlobM * cSys);
    problemData.cDiscRK{end, 2} = reshape(cSys(:, 1), [], problemData.g.numT).';
    problemData.cDiscRK{end, 3} = reshape(cSys(:, 2), [], problemData.g.numT).';
  end % if
  
  % Determine smoothed height
  [Q,~] = quadRule1D(problemData.qOrd);
  problemData.hSmoothV0T1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,end), [4 3], 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), [1 2], 2);
  assert(max(max(abs(hSmoothV0T1D - problemData.hSmoothV0T1D))) < 1e-12, 'hSmoothV0T1D and problemData.heightV0T1D must be the same');
  problemData.hSmoothQ0T1D = problemData.hSmoothV0T1D(:,1) * (1-Q) + problemData.hSmoothV0T1D(:,2) * Q;
end % if
end % function

 
% 
%> @brief Assembles static matrices
%>
%> This routine is called by adaptFreeSurface() to (re-)assemble static matrices
%> after modifying the mesh.
%>
%> The following matrices are assembled:
%> - @f$\mathsf{M}@f$
%> - @f$\check{\mathsf{H}}@f$
%> - @f$\check{\mathsf{Q}}@f$
%> - @f$\check{\mathsf{Q}}_\mathrm{bdr}@f$
%> - @f$\mathsf{H}@f$
%> - @f$\mathsf{Q}@f$
%> - @f$\mathsf{Q}_\mathrm{bdr}@f$
%> - @f$\mathsf{Q}_\mathrm{top}@f$
%> - @f$\mathsf{Q}_\mathrm{avg}@f$
%> - @f$\mathsf{Q}_\mathrm{up}@f$
%> - @f$\mathsf{S}_u@f$
%> - @f$\mathsf{S}_h@f$
%> - @f$\bar{\mathsf{S}}_h@f$
%>
%> For a definition of these matrices, see @ref RRAFK2018 .
%>
%> @param  problemData  A struct with problem parameters, basis functions and
%>                      representation matrices of the current solution
%>                      @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description. @f$[\text{struct}]@f$
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
function problemData = assembleStaticMatrices(problemData)
% Mass matrix (I, VIII)
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);

%% Momentum equation
% Element integral for height-term in momentum equation (III)
globVeeH = assembleMatElemTetraDphiPhi1D(problemData.g, problemData.g.g1D, problemData.hatVeeH);

% Interior edge integral for height-term in momentum equation (VI)
markE0T = problemData.g.markE0Tint | (problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrH);
globVeeQ = assembleMatEdgeTetraPhiPhi1DNu(problemData.g, problemData.g.g1D, markE0T, problemData.hatVeeQdiag, problemData.hatVeeQoffdiag);

% Boundary edge integrals
markE0T = problemData.g.markE0TbdrTop | problemData.g.markE0TbdrBot | (problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrH);
globVeeQbdr = assembleMatEdgeTetraPhiIntPhi1DIntNu(problemData.g, problemData.g.g1D, markE0T, problemData.hatVeeQdiag);

% Combine matrices
problemData.globVeeHQ = problemData.gConst * (globVeeH{1} - globVeeQ{1} - globVeeQbdr{1});

%% Flux equation
% Element integral in flux and continuity equation (IX, XI)
globH = assembleMatElemDphiPhi(problemData.g, problemData.hatH);

% Boundary edge integral without Dirichlet data for U in flux and continuity equation (X, XII)
markE0T = problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrU;
globQbdr = assembleMatEdgePhiIntPhiIntNu(problemData.g, markE0T, problemData.hatQdiag);
globQtop = assembleMatEdgePhiIntPhiIntNu(problemData.g, problemData.g.markE0TbdrTop, problemData.hatQdiag);

% Interior edge integral in flux equation (X)
globQ = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint, problemData.hatQdiag, problemData.hatQoffdiag);

% Combine matrices
problemData.globHQ = cellfun(@(H, Q, Qbdr, Qtop) H - Q - Qbdr - Qtop, globH, globQ, globQbdr, globQtop, 'UniformOutput', false);

%% Continuity equation
% Horizontal interior edge integral with first normal component in continuity equation (XII)
globQAvg = assembleMatEdgePhiPhiNu(problemData.g, problemData.g.markE0Tint & problemData.g.markE0Th, problemData.hatQdiag, problemData.hatQoffdiag);

% Horizontal interior and top boundary edge integral with second normal component in continuity equation (XII)
markE0T = (problemData.g.markE0Tint & problemData.g.markE0Th) | problemData.g.markE0TbdrTop;
globQup = problemData.fn_assembleMatEdgeTetraHorizPhiPhiNuBottomUp(problemData.g, markE0T, problemData.hatQdiag, problemData.hatQoffdiag);

% Combine matrices
problemData.globHQup = globH{2} - globQup;
problemData.globHQavg = -globH{1} + globQAvg{1} + globQtop{1};

%% Helper matrices for assembly of jump terms in Lax-Friedrichs Riemann solver
% Helper matrix for jumps over vertical edges in momentum equation
markE0T = (problemData.g.markE0Tint & problemData.g.markE0Tv) | (problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU);
problemData.globSu = assembleMatEdgePhiIntPerQuad(problemData.g, markE0T, problemData.hatSdiag);

% Helper matrix for jumps over vertical edges in continuity equation
markE0T = (problemData.g.markE0Tint & problemData.g.markE0Tv) | (problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrH);
problemData.globSh = assembleMatEdgePhiIntPerQuad(problemData.g, markE0T, problemData.hatSdiag);

% Helper matrix for jumps over vertices in free surface equation 
problemData.globBarS = assembleMatEdgeTetraPhi1DIntPerQuad(problemData.g, problemData.g.g1D, markE0T, problemData.hatBarSdiag);
end % function