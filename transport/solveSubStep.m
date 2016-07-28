function problemData = solveSubStep(problemData, ~, nSubStep)
K = problemData.K;
N = problemData.N;
p = problemData.p;

qOrd2D = max(2* p,1); 
[Q1, Q2, W] = quadRule2D(qOrd2D); 
numQuad2D = length(W);

% Assembly of time-dependent global matrices
globG = assembleMatElemDphiPhiFuncDiscVec(problemData.g, problemData.hatG, problemData.u1Disc, problemData.u2Disc);
globR = assembleMatEdgePhiPhiValUpwind(problemData.g, problemData.g.markE0Tint, problemData.hatRdiagOnQuad, ...
                                       problemData.hatRoffdiagOnQuad,  problemData.vNormalOnQuadEdge, problemData.g.areaE0TbdrNotN);
% Building the system
sysA = -globG{1} - globG{2} + globR;

for species = 1:problemData.numSpecies
  % L2 projections of algebraic coefficients
  fDisc  = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fCont{species}(problemData.timeLvls(nSubStep),x1,x2), ...
                                    2*p, problemData.hatM, problemData.basesOnQuad);

  % Assembly of Dirichlet boundary contributions
  globKD = assembleVecEdgePhiIntFuncContVal(problemData.g, problemData.g.markE0TbdrD, ...
            @(x1,x2) problemData.cDCont{species}(problemData.timeLvls(nSubStep),x1,x2),  problemData.vNormalOnQuadEdge, N, ...
            problemData.basesOnQuad, problemData.g.areaE0TbdrD);

  % Assembly of Neumann boundary contributions
  gNUpwind = @(x1,x2) (problemData.gNCont{species}(problemData.timeLvls(nSubStep),x1,x2) <= 0) .* problemData.gNCont{species}(problemData.timeLvls(nSubStep),x1,x2);
  globKN = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
            gNUpwind, N, problemData.basesOnQuad);

  % Assembly of the source contribution
  globL = problemData.globM * reshape(fDisc', K*N, 1) ...
        + problemData.globT * reshape(problemData.reactions{species}(problemData.timeLvls(nSubStep), problemData.g.mapRef2Phy(1, Q1, Q2), problemData.g.mapRef2Phy(2, Q1, Q2), problemData.cQ0T).', numQuad2D*K, 1);

  % right hand side
  sysV = globL - globKD - globKN;

  % Computing the discrete time derivative
  cDiscDot = problemData.globM \ (sysV - sysA * problemData.cDiscRK{species});

  % Apply slope limiting to time derivative
  if problemData.isSlopeLim{species}
    cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot, [N K])', problemData.globM, problemData.globMDiscTaylor);
    cDiscDotTaylorLim = applySlopeLimiterTaylor(problemData.g, cDiscDotTaylor, problemData.g.markV0TbdrD, NaN(K,3), problemData.basesOnQuad, problemData.typeSlopeLim{species});
    cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + problemData.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
    cDiscDot = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', problemData.globM, problemData.globMDiscTaylor)', [K*N 1]);
  end % if

  % Compute next step
  problemData.cDiscRK{species} = problemData.omega(nSubStep) * problemData.cDiscRK0{species} + (1 - problemData.omega(nSubStep)) * (problemData.cDiscRK{species} + problemData.tau * cDiscDot);

  % Limiting the solution
  if problemData.isSlopeLim{species}
    cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont{species}(problemData.timeLvls(nSubStep), x1, x2));
    problemData.cDiscRK{species} = reshape(applySlopeLimiterDisc(problemData.g, reshape(problemData.cDiscRK{species}, [N K])', problemData.g.markV0TbdrD, ...
                                   cDV0T, problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim{species})', [K*N 1]);
  end % if
end % for

end % function
