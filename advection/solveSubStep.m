function problemData = solveSubStep(problemData, ~, nSubStep)
K = problemData.K;
N = problemData.N;

% L2 projections of algebraic coefficients
fDisc  = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fCont(problemData.t(nSubStep),x1,x2),  ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);
u1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u1Cont(problemData.t(nSubStep),x1,x2), ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);
u2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u2Cont(problemData.t(nSubStep),x1,x2), ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);

% Evaluate normal velocity in quadrature points of edges
vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, @(x1,x2) problemData.u1Cont(problemData.t(nSubStep),x1,x2), ...
                      @(x1,x2) problemData.u2Cont(problemData.t(nSubStep),x1,x2), 2*problemData.p+1);

% Assembly of time-dependent global matrices
globG = assembleMatElemDphiPhiFuncDiscVec(problemData.g, problemData.hatG, u1Disc, u2Disc);
globR = assembleMatEdgePhiPhiValUpwind(problemData.g, problemData.hatRdiagOnQuad, ...
                                       problemData.hatRoffdiagOnQuad, vNormalOnQuadEdge);

% Assembly of Dirichlet boundary contributions
globKD = assembleVecEdgePhiIntFuncContVal(problemData.g, problemData.g.markE0TbdrD, ...
          @(x1,x2) problemData.cDCont(problemData.t(nSubStep),x1,x2), vNormalOnQuadEdge, N, ...
          problemData.basesOnQuad, problemData.g.areaE0TbdrD);

% Assembly of Neumann boundary contributions
gNUpwind = @(x1,x2) (problemData.gNCont(problemData.t(nSubStep),x1,x2) <= 0) .* problemData.gNCont(problemData.t(nSubStep),x1,x2);
globKN = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
          gNUpwind, N, problemData.basesOnQuad);

% Assembly of the source contribution
globL = problemData.globM * reshape(fDisc', K*N, 1);

% Building the system
sysA = -globG{1} - globG{2} + globR;
sysV = globL - globKD - globKN;

% Computing the discrete time derivative
cDiscDot = problemData.globM \ (sysV - sysA * problemData.cDiscRK{nSubStep});

% Apply slope limiting to time derivative
if problemData.isSlopeLim
  cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot, [N K])', problemData.globM, problemData.globMDiscTaylor);
  cDiscDotTaylorLim = applySlopeLimiterTaylor(problemData.g, cDiscDotTaylor, problemData.g.markV0TbdrD, NaN(K,3), problemData.basesOnQuad, problemData.typeSlopeLim);
  cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + problemData.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
  cDiscDot = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', problemData.globM, problemData.globMDiscTaylor)', [K*N 1]);
end % if

% Compute next step
problemData.cDiscRK{nSubStep + 1} = problemData.omega(nSubStep) * problemData.cDiscRK{1} + (1 - problemData.omega(nSubStep)) * (problemData.cDiscRK{nSubStep} + problemData.tau * cDiscDot);

% Limiting the solution
if problemData.isSlopeLim
  cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont(problemData.t(nSubStep), x1, x2));
  problemData.cDiscRK{nSubStep + 1} = reshape(applySlopeLimiterDisc(problemData.g, reshape(problemData.cDiscRK{nSubStep + 1}, [N K])', problemData.g.markV0TbdrD, ...
                          cDV0T, problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim)', [K*N 1]);
end % if

end % function
