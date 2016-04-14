function problemData = outputStep(problemData, nStep)
K = problemData.K;
N = problemData.N;
%% Visualization
if problemData.isVisSol
  cDisc = reshape(problemData.sysY(2*K*N+1 : 3*K*N), N, K)';
  cLagr = projectDataDisc2DataLagr(cDisc);
  visualizeDataLagr(problemData.g, cLagr, 'c_h', problemData.outputBasename, nStep, problemData.outputTypes);
end % if
end % function

