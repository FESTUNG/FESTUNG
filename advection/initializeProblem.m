function problemData = initializeProblem(problemData)
%% Initial data.
problemData.cDisc = projectFuncCont2DataDisc(problemData.g, problemData.c0Cont, 2*problemData.p+1, ...
                                             problemData.hatM, problemData.basesOnQuad);
if problemData.isSlopeLim
  cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont(0, x1, x2));
  problemData.cDisc = applySlopeLimiterDisc(problemData.g, problemData.cDisc, ...
                        problemData.g.markV0TbdrD, cDV0T, problemData.globM, ...
                        problemData.globMDiscTaylor, problemData.basesOnQuad, ...
                        problemData.typeSlopeLim);
end % if
fprintf('L2 error w.r.t. the initial condition: %g\n', ...
  computeL2Error(problemData.g, problemData.cDisc, problemData.c0Cont, 2*problemData.p, problemData.basesOnQuad));
%% visualization of inital condition.
if problemData.isVisSol
  cLagrange = projectDataDisc2DataLagr(problemData.cDisc);
  visualizeDataLagr(problemData.g, cLagrange, 'u_h', problemData.outputBasename, 0, problemData.outputTypes)
end
fprintf('Starting time integration from 0 to %g using time step size %g (%d steps).\n', ...
  problemData.tEnd, problemData.tau, problemData.numSteps)
end % function