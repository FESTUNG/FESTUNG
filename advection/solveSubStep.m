function problemData = solveSubStep(problemData, ~, nSubStep)
K = problemData.K;
N = problemData.N;

% Building the system
sysA = -problemData.globG{1} - problemData.globG{2} + problemData.globR;
sysV = problemData.globL - problemData.globKD - problemData.globKN;

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
