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
  
  % Re-assemble static matrices
  problemData = assembleStaticMatrices(problemData);
  
  % Determine smoothed height
  [Q,~] = quadRule1D(problemData.qOrd);
  problemData.hSmoothV0T1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,end), [4 3], 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), [1 2], 2);
  assert(max(max(abs(hSmoothV0T1D - problemData.hSmoothV0T1D))) < 1e-12, 'hSmootV0T1D and problemData.heightV0T1D must be the same');
  problemData.hSmoothQ0T1D = problemData.hSmoothV0T1D(:,1) * (1-Q) + problemData.hSmoothV0T1D(:,2) * Q;
end % if
end % function
%
function problemData = assembleStaticMatrices(problemData)
% Mass matrix (I, VIII)
problemData.globM = assembleMatElemPhiPhi(problemData.g, problemData.hatM);

%% Momentum equation
% Element integral for height-term in momentum equation (III)
tildeGlobH = assembleMatElemTetraDphiPhi1D(problemData.g, problemData.g.g1D, problemData.tildeHatH);

% Interior edge integral for height-term in momentum equation (VI)
tildeGlobQ = assembleMatEdgeTetraPhiPhi1DNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tint, problemData.tildeHatQdiag, problemData.tildeHatQoffdiag);

% Boundary edge integral for height term in momentum equation with no prescribed Dirichlet-data (VI)
tildeGlobQbdr = assembleMatEdgeTetraPhiIntPhi1DIntNu(problemData.g, problemData.g.g1D, problemData.g.markE0Tbdr & ~(problemData.g.markE0TbdrH | problemData.g.markE0TbdrRiemH), problemData.tildeHatQdiag);

% Boundary edge integral for height term in momentum equation with Riemann solver and prescribed Dirichlet-data (VI)
tildeGlobQRiem = assembleMatEdgeTetraPhiIntPhi1DIntNu(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrRiemH, problemData.tildeHatQdiag);

% Combine matrices
problemData.tildeGlobHQ = problemData.gConst * (tildeGlobH{1} - tildeGlobQ{1} - tildeGlobQbdr{1} - 0.5 * tildeGlobQRiem{1});

%% Flux and continuity equation
% Element integral in flux and continuity equation (IX, XI)
globH = assembleMatElemDphiPhi(problemData.g, problemData.hatH);

% Boundary edge integral without Dirichlet data for U in flux and continuity equation (X, XII)
globQbdr = assembleMatEdgeTetraPhiIntPhiIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrU, problemData.hatQdiag);

%% Flux equation
% Interior edge integral in flux equation (X)
globQ = assembleMatEdgeTetraPhiPhiNu(problemData.g, problemData.g.markE0Tint, problemData.hatQdiag, problemData.hatQoffdiag);

% Combine matrices
problemData.globHQ = cellfun(@(H, Q, Qbdr) H - Q - Qbdr, globH, globQ, globQbdr, 'UniformOutput', false);

%% Continuity equation
% Horizontal interior edge integral with first normal component in continuity equation (XII)
globQAvg = assembleMatEdgeTetraPhiPhiNu(problemData.g, problemData.g.markE0Tint & problemData.g.markE0Th, problemData.hatQdiag, problemData.hatQoffdiag);

% Horizontal interior and top boundary edge integral with second normal component in continuity equation (XII)
globQup = problemData.fn_assembleMatEdgeTetraHorizPhiPhiNuBottomUp(problemData.g, (problemData.g.markE0Tint | problemData.g.markE0TbdrTop) & problemData.g.markE0Th, problemData.hatQdiag, problemData.hatQoffdiag);

% Combine matrices
problemData.globHQup = globH{2} - globQup - globQAvg{2};
problemData.globHQavg = -globH{1} + globQAvg{1} + globQbdr{1};

%% Helper matrices for assembly of jump terms in Lax-Friedrichs Riemann solver
% Helper matrix for jumps over vertical interior edges in momentum and continuity equation (VI, XII)
problemData.globS = assembleMatEdgeTetraPhiIntPerQuad(problemData.g, problemData.g.markE0Tint & problemData.g.markE0Tv, problemData.hatSdiag);

% Helper matrix for jumps over vertical boundary edges in momentum equation with prescribed Dirichlet data for u and Riemann solver (VI)
problemData.globSuRiem = assembleMatEdgeTetraPhiIntPerQuad(problemData.g, problemData.g.markE0TbdrRiemU, problemData.hatSdiag);

% Helper matrix for jumps over vertical boundary edges in continuity equation with prescribed Dirichlet data for h and Riemann solver (VI)
problemData.globShRiem = assembleMatEdgeTetraPhiIntPerQuad(problemData.g, problemData.g.markE0TbdrRiemH, problemData.hatSdiag);

% Helper matrix for jumps over interior vertices in free surface equation (XV)
problemData.barGlobS = assembleMatEdgeTetraPhi1DIntPerQuad(problemData.g, problemData.g.g1D, problemData.g.markE0Tint, problemData.barHatSdiag);

% Helper matrix for jumps over boundary vertices in free surface equation with prescribed Dirichlet data for h and Riemann solver (XV)
problemData.barGlobSRiem = assembleMatEdgeTetraPhi1DIntPerQuad(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrRiemH, problemData.barHatSdiag);
end % function