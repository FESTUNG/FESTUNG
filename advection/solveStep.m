function problemData = solveStep(problemData, nStep)
K = problemData.K;
N = problemData.N;

% Obtain Runge-Kutta rule
[t, omega] = rungeKuttaSSP(problemData.ordRK, problemData.tau, (nStep - 1) * problemData.tau);
cDiscRK = cell(length(omega)+1, 1); 
cDiscRK{1} = reshape(problemData.cDisc', [K*N 1]);

%% Perform Runge-Kutta steps
for rkStep = 1 : length(omega)
  % L2 projections of algebraic coefficients
  fDisc  = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fCont(t(rkStep),x1,x2),  ...
                                    2*problemData.p, problemData.hatM, problemData.basesOnQuad);
  u1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u1Cont(t(rkStep),x1,x2), ...
                                    2*problemData.p, problemData.hatM, problemData.basesOnQuad);
  u2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u2Cont(t(rkStep),x1,x2), ...
                                    2*problemData.p, problemData.hatM, problemData.basesOnQuad);
                                  
  % Evaluate normal velocity in quadrature points of edges
  vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, @(x1,x2) problemData.u1Cont(t(rkStep),x1,x2), ...
                        @(x1,x2) problemData.u2Cont(t(rkStep),x1,x2), 2*problemData.p+1);
                      
  % Assembly of time-dependent global matrices
  globG = assembleMatElemDphiPhiFuncDiscVec(problemData.g, problemData.hatG, u1Disc, u2Disc);
  globR = assembleMatEdgePhiPhiValUpwind(problemData.g, ~problemData.g.markE0TbdrN, ...
                                         problemData.hatRdiagOnQuad, problemData.hatRoffdiagOnQuad, ...
                                         vNormalOnQuadEdge, problemData.g.areaE0TbdrNotN);

  % Assembly of Dirichlet boundary contributions
  globKD = assembleVecEdgePhiIntFuncContVal(problemData.g, problemData.g.markE0TbdrD, ...
            @(x1,x2) problemData.cDCont(t(rkStep),x1,x2), vNormalOnQuadEdge, N, ...
            problemData.basesOnQuad, problemData.g.areaE0TbdrD);
          
  % Assembly of Neumann boundary contributions
  gNUpwind = @(x1,x2) (problemData.gNCont(t(rkStep),x1,x2) <= 0) .* problemData.gNCont(t(rkStep),x1,x2);
  globKN = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
            gNUpwind, N, problemData.basesOnQuad);
          
  % Assembly of the source contribution
  globL = problemData.globM * reshape(fDisc', K*N, 1);
  
  % Building the system
  sysA = -globG{1} - globG{2} + globR;
  sysV = globL - globKD - globKN;
  
  % Computing the discrete time derivative
  cDiscDot = problemData.globM \ (sysV - sysA * cDiscRK{rkStep});
  
  % Apply slope limiting to time derivative
  if problemData.isSlopeLim
    cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot, [N K])', problemData.globM, problemData.globMDiscTaylor);
    cDiscDotTaylorLim = applySlopeLimiterTaylor(problemData.g, cDiscDotTaylor, problemData.g.markV0TbdrD, NaN(K,3), problemData.basesOnQuad, problemData.typeSlopeLim);
    cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + problemData.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
    cDiscDot = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', problemData.globM, problemData.globMDiscTaylor)', [K*N 1]);
  end
  
  % Compute next step
  cDiscRK{rkStep + 1} = omega(rkStep) * cDiscRK{1} + (1 - omega(rkStep)) * (cDiscRK{rkStep} + problemData.tau * cDiscDot);
  
  % Limiting the solution
  if problemData.isSlopeLim
    cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont(t(rkStep), x1, x2));
    cDiscRK{rkStep + 1} = reshape(applySlopeLimiterDisc(problemData.g, reshape(cDiscRK{rkStep + 1}, [N K])', problemData.g.markV0TbdrD, ...
                            cDV0T, problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim)', [K*N 1]);
  end % if
end % for

%% Reshape and store solution in problemData
problemData.cDisc = reshape(cDiscRK{end}, N, K)';
end % function