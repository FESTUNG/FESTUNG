function problemData = initializeProblem(problemData)

problemData.isFinished = false;

%% Initial data.
h0Cont = @(x,y) problemData.hSol(0,x,y);
hDisc = projectFuncCont2DataDisc(problemData.g, h0Cont, ...
                      problemData.p, problemData.ord, problemData.hatMc, problemData.basesOnQuad);
problemData.sysY = [ zeros(2 * problemData.NT * problemData.N, 1) ; ...
                     reshape(hDisc', problemData.NT * problemData.N, 1) ];
fprintf('L2 error w.r.t. the initial condition: %g\n', ...
  computeL2Error(problemData.g, hDisc, h0Cont, problemData.ord, problemData.basesOnQuad));

%% visualization of inital condition.
if problemData.isVisSol
  hLagrange = projectDataDisc2DataLagr(hDisc);
  visualizeDataLagr(problemData.g, hLagrange, 'h_h', problemData.outputBasename, 0);
end % if

end % function