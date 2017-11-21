function problemData = initializeProblem(problemData)
problemData.isFinished = false;
K = problemData.g.numT;
N = problemData.N;

%% Initial data.
h0Cont = @(x1,x2) problemData.hCont(0,x1,x2);
q10Cont = @(x1,x2) problemData.q1Cont(0,x1,x2);
q20Cont = @(x1,x2) problemData.q2Cont(0,x1,x2);

hDisc = projectFuncCont2DataDiscTetra(problemData.g, h0Cont, problemData.qOrd, ...
                                      problemData.globM, problemData.basesOnQuad);
q1Disc = projectFuncCont2DataDiscTetra(problemData.g, q10Cont, problemData.qOrd, ...
                                       problemData.globM, problemData.basesOnQuad);
q2Disc = projectFuncCont2DataDiscTetra(problemData.g, q20Cont, problemData.qOrd, ...
                                       problemData.globM, problemData.basesOnQuad);

problemData.sysY = [ reshape(q1Disc', K * N, 1) ; ...
                     reshape(q2Disc', K * N, 1) ; ...
                     reshape(hDisc', K * N, 1) ];

%% Error computation and visualization of inital condition.
fprintf('L2 error w.r.t. the initial condition: %g\n', ...
  computeL2ErrorTetra(problemData.g, hDisc, h0Cont, problemData.qOrd+1, problemData.basesOnQuad));

if problemData.isVisSol
  cLagr = { projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(2*K*N+1 : 3*K*N), N, K)'), ...
            projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(1 : K*N), N, K)'), ...
            projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(K*N+1 : 2*K*N), N, K)') };
  visualizeDataLagrTetra(problemData.g, cLagr, {'h', 'q1', 'q2'}, problemData.outputBasename, 0, problemData.outputTypes, struct('q', {{'q1','q2'}}));
end % if

end % function