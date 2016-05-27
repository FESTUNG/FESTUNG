function problemData = initializeProblem(problemData)
problemData.isFinished = false;
%% Initial data.
cDisc = projectFuncCont2DataDisc(problemData.g, problemData.c0Cont, ...
                      2 * problemData.p, problemData.hatM, problemData.basesOnQuad);
problemData.sysY = [ zeros(2 * problemData.K * problemData.N, 1) ; ...
                     reshape(cDisc', problemData.K * problemData.N, 1) ];
fprintf('L2 error w.r.t. the initial condition: %g\n', ...
  computeL2Error(problemData.g, cDisc, problemData.c0Cont, 2 * problemData.p, problemData.basesOnQuad));
%% visualization of inital condition.
if problemData.isVisSol
  cLagrange = projectDataDisc2DataLagr(cDisc);
  visualizeDataLagr(problemData.g, cLagrange, 'u_h', problemData.outputBasename, 0, problemData.outputTypes)
end % if
end % function