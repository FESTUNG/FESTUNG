function problemData = initializeProblem(problemData)
problemData.isFinished = false;

%% Initial data.
h0Cont = @(x1,x2) problemData.hCont(0,x1,x2);
q10Cont = @(x1,x2) problemData.q1Cont(0,x1,x2);
q20Cont = @(x1,x2) problemData.q2Cont(0,x1,x2);

hDisc = projectFuncCont2DataDiscTetra(problemData.g, h0Cont, problemData.N, problemData.qOrd, ...
                                      problemData.globM, problemData.basesOnQuad);
q1Disc = projectFuncCont2DataDiscTetra(problemData.g, q10Cont, problemData.N, problemData.qOrd, ...
                                       problemData.globM, problemData.basesOnQuad);
q2Disc = projectFuncCont2DataDiscTetra(problemData.g, q20Cont, problemData.N, problemData.qOrd, ...
                                       problemData.globM, problemData.basesOnQuad);

problemData.sysY = [ reshape(q1Disc', problemData.g.numT * problemData.N, 1) ; ...
                     reshape(q2Disc', problemData.g.numT * problemData.N, 1) ; ...
                     reshape(hDisc', problemData.g.numT * problemData.N, 1) ];

%% Error computation and visualization of inital condition.
fprintf('L2 error w.r.t. the initial condition: %g\n', ...
  computeL2ErrorTetra(problemData.g, hDisc, h0Cont, problemData.qOrd+1, problemData.basesOnQuad));

if problemData.isVisSol
  hLagr = projectDataDisc2DataLagrTensorProduct(hDisc);
  visualizeDataLagrTetra(problemData.g, hLagr, 'h', problemData.outputBasename, 0);
end % if

end % function